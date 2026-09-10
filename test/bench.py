#!/usr/bin/env python3
"""Before/after benchmark: adze vs. a reference ADZE 1.0 build.

Measures the axes the code review identified as the scaling problems -- locus
filtering, the MAX_G sweep, the tuple enumeration, and the read pass's time and
memory -- plus thread scaling of the reworked statistics passes.

    ./test/bench.py --candidate src/adze --reference /path/to/adze-1.0 \
        --out bench.json

Every timing is the best of --reps runs (default 2); memory is the peak
resident set size reported by wait4.
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen_data  # noqa: E402


class Runner:
    def __init__(self, workdir, reps):
        self.workdir = workdir
        self.reps = reps
        os.makedirs(workdir, exist_ok=True)

    def time(self, binary, paramfile, extra=()):
        best = float("inf")
        for _ in range(self.reps):
            t0 = time.perf_counter()
            proc = subprocess.run([os.path.abspath(binary), paramfile, *extra],
                                  cwd=self.workdir, capture_output=True, text=True)
            assert proc.returncode == 0, (binary, paramfile, proc.stderr[-400:],
                                          proc.stdout[-400:])
            best = min(best, time.perf_counter() - t0)
        return best

    def peak_rss(self, binary, paramfile, extra=()):
        peak = 0
        for _ in range(self.reps):
            proc = subprocess.Popen([os.path.abspath(binary), paramfile, *extra],
                                    cwd=self.workdir,
                                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            _, _, ru = os.wait4(proc.pid, 0)
            # bytes on macOS, KiB on Linux
            rss = ru.ru_maxrss if sys.platform == "darwin" else ru.ru_maxrss * 1024
            peak = max(peak, rss)
        return peak


def filtering(run, cand, ref):
    """Locus-filtering phase, isolated by differencing a no-op tolerance."""
    rows = []
    for nloci in [1000, 2000, 4000, 8000]:
        d = os.path.join(run.workdir, "filt_%d.stru" % nloci)
        meta = gen_data.write_dataset(d, 10, 25, nloci, 10, missing=0.30, seed=11)
        keep = os.path.join(run.workdir, "filt_keep.par")
        drop = os.path.join(run.workdir, "filt_drop.par")
        gen_data.write_paramfile(keep, meta, os.path.basename(d), "fk", g=2, tol=0.9)
        gen_data.write_paramfile(drop, meta, os.path.basename(d), "fd", g=2, tol=0.38)

        row = {"loci": nloci}
        for tag, binary in (("ref", ref), ("cand", cand)):
            base = run.time(binary, "filt_keep.par", ("-tnocalc", "1"))
            full = run.time(binary, "filt_drop.par", ("-tnocalc", "1"))
            row[tag] = full - base
        with open(os.path.join(run.workdir, "fd_p_deletedloci")) as fh:
            row["deleted"] = int(fh.readline().split()[0])
        rows.append(row)
        print("  filtering loci=%-5d ref %6.2fs cand %6.2fs" %
              (nloci, row["ref"], row["cand"]))
    return rows


def maxg(run, cand, ref):
    """Statistics phase vs. MAX_G, isolated by differencing the MAX_G 2 run."""
    d = os.path.join(run.workdir, "maxg.stru")
    meta = gen_data.write_dataset(d, 10, 50, 1000, 20, seed=3)
    rows = []
    base = {}
    for G in [2, 10, 20, 40, 80]:
        par = os.path.join(run.workdir, "maxg_%d.par" % G)
        gen_data.write_paramfile(par, meta, "maxg.stru", "mg", g=G, tol=1)
        row = {"max_g": G}
        for tag, binary in (("ref", ref), ("cand", cand)):
            t = run.time(binary, os.path.basename(par))
            if G == 2:
                base[tag] = t
            row[tag + "_total"] = t
            row[tag] = t - base[tag]
        rows.append(row)
        print("  MAX_G %-4d ref %6.2fs cand %6.2fs (statistics phase)" %
              (G, row["ref"], row["cand"]))
    return rows


def tuples(run, cand, ref):
    """Tuple phase vs. the number of tuples evaluated."""
    d = os.path.join(run.workdir, "tup.stru")
    meta = gen_data.write_dataset(d, 12, 10, 500, 15, seed=5)
    par0 = os.path.join(run.workdir, "tup_base.par")
    gen_data.write_paramfile(par0, meta, "tup.stru", "tb", g=10, tol=1)
    base = {tag: run.time(b, "tup_base.par") for tag, b in (("ref", ref), ("cand", cand))}

    rows = []
    for k in [1, 2, 3, 4]:
        par = os.path.join(run.workdir, "tup_%d.par" % k)
        gen_data.write_paramfile(par, meta, "tup.stru", "tk", g=10, tol=1,
                                 comb=1, k=str(k))
        row = {"k": k, "tuples": math.comb(12, k)}
        for tag, binary in (("ref", ref), ("cand", cand)):
            row[tag] = run.time(binary, os.path.basename(par)) - base[tag]
        rows.append(row)
        print("  k=%d (%4d tuples) ref %6.2fs cand %6.2fs" %
              (k, row["tuples"], row["ref"], row["cand"]))
    return rows


def reading(run, cand, ref):
    """Read pass: wall time and peak memory vs. input size."""
    rows = []
    for nloci, nind in [(2000, 25), (4000, 25), (8000, 25), (8000, 100)]:
        d = os.path.join(run.workdir, "read_%d_%d.stru" % (nloci, nind))
        meta = gen_data.write_dataset(d, 10, nind, nloci, 10, seed=13)
        par = os.path.join(run.workdir, "read.par")
        gen_data.write_paramfile(par, meta, os.path.basename(d), "rd", g=2, tol=1)
        row = {"loci": nloci, "rows": meta["dlines"],
               "cells": meta["dlines"] * nloci,
               "file_mb": os.path.getsize(d) / 1e6}
        for tag, binary in (("ref", ref), ("cand", cand)):
            row[tag + "_time"] = run.time(binary, "read.par")
            row[tag + "_rss"] = run.peak_rss(binary, "read.par")
        rows.append(row)
        print("  read %5d loci x %4d rows (%5.1f MB): ref %5.2fs/%6.1f MiB  "
              "cand %5.2fs/%6.1f MiB" %
              (nloci, row["rows"], row["file_mb"], row["ref_time"],
               row["ref_rss"] / 2**20, row["cand_time"], row["cand_rss"] / 2**20))
    return rows


def threads(run, cand, max_threads=8):
    """Thread scaling of the statistics phase (candidate only)."""
    d = os.path.join(run.workdir, "thr.stru")
    meta = gen_data.write_dataset(d, 10, 100, 2000, 20, seed=3)
    read_par = os.path.join(run.workdir, "thr_read.par")
    work_par = os.path.join(run.workdir, "thr_work.par")
    gen_data.write_paramfile(read_par, meta, "thr.stru", "tr", g=2, tol=1)
    gen_data.write_paramfile(work_par, meta, "thr.stru", "tw", g=150, tol=1)

    serial_read = run.time(cand, "thr_read.par", ("--threads", "1"))
    rows = []
    n = 1
    while n <= max_threads:
        total = run.time(cand, "thr_work.par", ("--threads", str(n)))
        digest = None
        for name in ("tw_r", "tw_p"):
            with open(os.path.join(run.workdir, name), "rb") as fh:
                data = fh.read()
            digest = (digest or b"") + data
        import hashlib
        rows.append({"threads": n, "total": total,
                     "stats": total - serial_read,
                     "md5": hashlib.md5(digest).hexdigest()})
        print("  threads %d: statistics phase %5.2fs  md5 %s" %
              (n, rows[-1]["stats"], rows[-1]["md5"][:12]))
        n *= 2
    return {"read_only": serial_read, "rows": rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--workdir", default=os.path.join(HERE, "benchwork"))
    ap.add_argument("--out", default="bench.json")
    ap.add_argument("--reps", type=int, default=2)
    ap.add_argument("--max-threads", type=int, default=8)
    args = ap.parse_args()

    shutil.rmtree(args.workdir, ignore_errors=True)
    run = Runner(args.workdir, args.reps)

    results = {"platform": sys.platform, "reps": args.reps}
    print("locus filtering:")
    results["filtering"] = filtering(run, args.candidate, args.reference)
    print("MAX_G sweep:")
    results["maxg"] = maxg(run, args.candidate, args.reference)
    print("tuple enumeration:")
    results["tuples"] = tuples(run, args.candidate, args.reference)
    print("reading:")
    results["reading"] = reading(run, args.candidate, args.reference)
    print("thread scaling:")
    results["threads"] = threads(run, args.candidate, args.max_threads)

    with open(args.out, "w") as fh:
        json.dump(results, fh, indent=1)
    print("\nwrote %s" % args.out)
    shutil.rmtree(args.workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
