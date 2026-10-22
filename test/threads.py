#!/usr/bin/env python3
"""The thread count must not change any output.

Every case is run at each thread count in THREADS, in its own scratch
directory, and every file the run writes must be byte-identical to the
single-threaded run. The one exception is the run summary: it records the
elapsed time of each step, and the THREADS value it was given, and those two
kinds of line are removed before it is compared.

The cases cover each input path (VCF, gzipped VCF, STRUCTURE through its
count conversion, phased VCF in haplotype mode) and each kind of output
(statistics, _fulldata, combinations, both window types, --at-g, --legacy,
the deleted-loci list and the _blocks file). Some runs also force the work
into small batches, so that a batch boundary falls inside a window, inside a
haplotype block run and between combination batches:

  ADZE_SWEEP_BATCH   loci per batch in the sweep
  ADZE_SCAN_BATCH    records per batch in the counting pass over a VCF
  ADZE_TUPLE_BUDGET  bytes for the combination accumulators, which forces
                     the combinations into several passes

Thread counts include 3, which divides nothing evenly, and 64, which is more
threads than some batches have loci.

Two controls keep the suite from passing vacuously: every case must write at
least one non-empty result file, and a run with a different MAX_G must be
reported as different by the same comparison.

    ./test/threads.py --candidate src/adze [--workdir DIR] [-v]

Exit status is 0 only if every case passes.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_data

THREADS = (1, 2, 3, 4, 8, 64)

# Lines of the run summary that legitimately depend on the run rather than
# the data: elapsed times, and the thread count itself.
VOLATILE = re.compile(r"\(d:h:m:s\)|^THREADS\b|^--threads\b")


class Suite(object):
    def __init__(self, verbose):
        self.verbose = verbose
        self.passed = 0
        self.failed = []

    def check(self, name, ok, detail=""):
        if ok:
            self.passed += 1
            if self.verbose:
                print("PASS  %s" % name)
        else:
            self.failed.append((name, detail))
            print("FAIL  %s: %s" % (name, detail))


def fixtures(w):
    """Write the input files once; return {name: argument list}."""
    os.makedirs(w, exist_ok=True)
    snp = os.path.join(w, "snp")
    if not os.path.exists(snp + ".vcf"):
        gen_data.write_pair(snp, nloci=300, npops=5, nind=[6, 8, 8, 10, 12],
                            nall=2, missing=0.05, seed=11, chroms=2,
                            spacing=1000)
    msat = os.path.join(w, "msat")
    if not os.path.exists(msat + ".vcf"):
        gen_data.write_pair(msat, nloci=120, npops=6, nind=[5, 7, 9, 9, 11, 14],
                            nall=12, missing=0.08, seed=12, chroms=2,
                            spacing=5000)
    hap = os.path.join(w, "hap")
    if not os.path.exists(hap + ".vcf"):
        gen_data.write_phased(hap, nloci=400, npops=4, nind=10, nall=2,
                              missing=0.02, seed=13, chroms=2, spacing=500)
    return {
        "vcf": ["--data", "snp.vcf", "--samples", "snp.samples"],
        "vcfgz": ["--data", "snp.vcf.gz", "--samples", "snp.samples"],
        "stru": ["--data", "msat.stru", "--group-col", "2",
                 "--loci-map", "msat.map"],
        "msatvcf": ["--data", "msat.vcf", "--samples", "msat.samples"],
        "hap": ["--data", "hap.vcf", "--samples", "hap.samples"],
    }


# (name, input, extra arguments, extra environment)
CASES = [
    ("vcf.plain", "vcf", [], {}),
    ("vcf.full", "vcf", ["--max-g", "6", "--full-richness", "--full-private"], {}),
    ("vcf.full.batch5", "vcf", ["--max-g", "6", "--full-richness",
                                "--full-private"], {"ADZE_SWEEP_BATCH": "5"}),
    ("vcf.tolerance", "vcf", ["--tolerance", "0.02"], {}),
    ("vcf.windows_bp", "vcf", ["--window-bp", "20000", "--step-bp", "5000",
                               "--max-g", "8"], {"ADZE_SWEEP_BATCH": "7"}),
    ("vcf.windows_loci", "vcf", ["--window-loci", "25", "--step-loci", "10",
                                 "--min-window-loci", "3"], {}),
    ("vcf.scan_batch", "vcf", ["--tolerance", "0.04", "--window-bp", "15000",
                               "--step-bp", "5000"],
     {"ADZE_SCAN_BATCH": "3", "ADZE_SWEEP_BATCH": "4"}),
    ("vcf.at_g_max", "vcf", ["--at-g", "max", "--full-richness"], {}),
    ("vcf.legacy", "vcf", ["--legacy", "--max-g", "5"], {}),
    ("vcfgz.pops", "vcfgz", ["--pops", "POP0,POP2,POP4", "--max-g", "6"], {}),
    ("stru.tuples", "stru", ["--combinations", "--tuples-k", "1-3",
                             "--full-tuples", "--max-g", "6"], {}),
    ("stru.tuples.passes", "stru", ["--combinations", "--tuples-k", "2-3",
                                    "--max-g", "6", "--window-loci", "20"],
     {"ADZE_TUPLE_BUDGET": "20000", "ADZE_SWEEP_BATCH": "3"}),
    ("stru.windows_full", "stru", ["--window-bp", "100000", "--step-bp",
                                   "50000", "--full-richness",
                                   "--full-private", "--max-g", "8"], {}),
    ("msatvcf.exclude", "msatvcf", ["--exclude-pops", "POP1",
                                    "--tolerance", "0.1"], {}),
    ("hap.snps", "hap", ["--hap-snps", "4", "--max-g", "6",
                         "--full-richness"], {}),
    ("hap.bp.windows", "hap", ["--hap-bp", "2000", "--window-bp", "20000",
                               "--step-bp", "10000", "--max-g", "6"],
     {"ADZE_SWEEP_BATCH": "2"}),
    ("hap.min_snps", "hap", ["--hap-snps", "6", "--min-hap-snps", "6",
                             "--combinations", "--tuples-k", "2"], {}),
]


def place(src, dst):
    """Make an input visible in a run's directory: a link where the system
    allows one, a copy where it does not (Windows, without the privilege)."""
    try:
        os.symlink(src, dst)
    except (OSError, NotImplementedError):
        shutil.copyfile(src, dst)


def run(binary, workdir, args, env):
    return subprocess.run([os.path.abspath(binary)] + args, cwd=workdir,
                          capture_output=True, text=True,
                          env=dict(os.environ, **env))


def results(workdir, inputs):
    """{file name: comparable bytes} for every file the run wrote."""
    out = {}
    for f in sorted(os.listdir(workdir)):
        if f in inputs:
            continue
        with open(os.path.join(workdir, f), "rb") as fh:
            raw = fh.read()
        if "summary" in f:
            text = raw.decode(errors="replace").splitlines()
            raw = "\n".join(l for l in text if not VOLATILE.search(l)).encode()
        out[f] = raw
    return out


def first_difference(a, b):
    la, lb = a.splitlines(), b.splitlines()
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            return "line %d: %r vs %r" % (i + 1, x[:60], y[:60])
    return "line %d: one file ends" % (min(len(la), len(lb)) + 1)


def one_case(s, adze, w, inputs, name, args, env):
    base = None
    for t in THREADS:
        d = os.path.join(w, "%s.t%d" % (name, t))
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        for f in inputs:
            place(os.path.join(w, f), os.path.join(d, f))
        r = run(adze, d, args + ["--quiet", "--out-prefix", "out",
                                 "--threads", str(t)], env)
        got = results(d, inputs)
        if t == THREADS[0]:
            base = (r.returncode, got)
            nonempty = [f for f, b in got.items() if b and "summary" not in f]
            s.check("%s.writes_results" % name, r.returncode == 0 and nonempty,
                    "exit %d, %d non-empty result files; %s" % (
                        r.returncode, len(nonempty), r.stderr.strip()[-200:]))
            continue
        problems = []
        if r.returncode != base[0]:
            problems.append("exit %d vs %d" % (r.returncode, base[0]))
        if set(got) != set(base[1]):
            problems.append("files differ: %s" % sorted(set(got) ^ set(base[1])))
        for f in sorted(set(got) & set(base[1])):
            if got[f] != base[1][f]:
                problems.append("%s %s" % (f, first_difference(base[1][f], got[f])))
        s.check("%s.threads%d" % (name, t), not problems, "; ".join(problems[:3]))
        if not problems:
            shutil.rmtree(d)
    shutil.rmtree(os.path.join(w, "%s.t%d" % (name, THREADS[0])), ignore_errors=True)


def control(s, adze, w, inputs):
    """The comparison must see a real difference when there is one."""
    runs = {}
    for g in ("5", "6"):
        d = os.path.join(w, "control.g%s" % g)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        for f in inputs:
            place(os.path.join(w, f), os.path.join(d, f))
        run(adze, d, inputs_args["vcf"] + ["--quiet", "--out-prefix", "out",
                                           "--max-g", g], {})
        runs[g] = results(d, inputs)
        shutil.rmtree(d)
    s.check("control.detects_difference",
            runs["5"].get("out.richness") != runs["6"].get("out.richness"),
            "MAX_G 5 and 6 compared equal")


inputs_args = {}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--workdir", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "work-threads"))
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    inputs_args.update(fixtures(args.workdir))
    inputs = set(f for f in os.listdir(args.workdir)
                 if os.path.isfile(os.path.join(args.workdir, f)))

    s = Suite(args.verbose)
    control(s, args.candidate, args.workdir, inputs)
    for name, src, extra, env in CASES:
        one_case(s, args.candidate, args.workdir, inputs, name,
                 inputs_args[src] + ["--tolerance", "1"] + extra
                 if "--tolerance" not in extra else inputs_args[src] + extra,
                 env)

    print("\n%d passed, %d failed  (%d cases x %d thread counts)"
          % (s.passed, len(s.failed), len(CASES), len(THREADS)))
    return 1 if s.failed else 0


if __name__ == "__main__":
    sys.exit(main())
