#!/usr/bin/env python3
"""Check the sliding-window statistics against what the layout implies.

There is no reference implementation to diff against, so every case here is an
invariant that must hold whatever the numbers are:

  whole            one window spanning a single-chromosome dataset reproduces
                   the genome-wide output exactly
  recompute        every window's mean, variance and standard error, computed
                   here from the per-locus *_fulldata columns, matches
  layout           windows tile as asked: locus membership, chromosome
                   boundaries, and each interior locus in exactly two windows
                   when the step is half the width
  units            evenly spaced loci give the same statistics whether the
                   window is expressed in basepairs or in loci
  minloci          --min-window-loci drops exactly the sparse windows
  formats          the same genotypes windowed from VCF and from
                   STRUCTURE-plus-map give identical output
  threads          window output does not depend on --threads
  refusals         unsorted positions, a chromosome in two blocks, a locus
                   with no coordinate, and STRUCTURE without a map are all
                   refused

    ./test/windows.py --candidate src/adze

Exit status is 0 only if every case passes.
"""
import argparse
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen_data                                          # noqa: E402

# Output is printed to six significant digits, so a mean recomputed from the
# per-locus columns can differ in the last one.
RTOL = 2e-5

# A variance recomputed from those same rounded columns is looser: it is built
# from deviations, which are small differences of similar numbers, so the
# rounding in each value is amplified by roughly value/deviation -- about 2e-5
# here. The check still catches a wrong locus set or a wrong divisor, both of
# which are off by percents, not parts per million.
VAR_RTOL = 1e-3


def run(binary, workdir, args, expect=0):
    r = subprocess.run([os.path.abspath(binary)] + args, cwd=workdir,
                       capture_output=True, text=True)
    return r


def read_windows(path):
    """Window rows as a list of dicts, using the header for column names."""
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    head = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        rows.append(dict(zip(head, line.split("\t"))))
    return rows


def read_fulldata(path):
    """Per-locus values as {(group, g): [floats]}, plus the locus names."""
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    head = lines[0].split()
    names = head[3:-3]
    out = {}
    for line in lines[1:]:
        f = line.split()
        n = int(f[2])
        out[(f[0], int(f[1]))] = [float(x) for x in f[3:3 + n]]
    return names, out


def close(a, b, rtol=RTOL):
    if a == b:
        return True
    scale = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / scale <= rtol


class Suite:
    def __init__(self, verbose):
        self.npass = 0
        self.failures = []
        self.verbose = verbose

    def check(self, label, ok, detail=""):
        if ok:
            self.npass += 1
            if self.verbose:
                print("PASS  %s" % label)
        else:
            self.failures.append((label, detail))
            print("FAIL  %s: %s" % (label, detail))


def case_whole(s, adze, w):
    """One window over a single-chromosome dataset == the genome-wide run."""
    gen_data.write_pair(os.path.join(w, "one"), npops=3, nind=6, nloci=20,
                        nall=5, missing=0.05, seed=11, chroms=1, spacing=1000)
    base = ["--data", "one.vcf", "--samples", "one.samples", "--max-g", "5"]
    run(adze, w, base + ["--out-prefix", "one_g"])
    run(adze, w, base + ["--out-prefix", "one_w", "--window-bp", "10000000"])

    genome = {}
    for line in open(os.path.join(w, "one_g.richness")):
        f = line.split()
        if f:
            genome[(f[0], int(f[1]))] = f[2:]

    rows = read_windows(os.path.join(w, "one_w.richness_windows"))
    s.check("whole.rowcount", len(rows) == len(genome),
            "%d window rows vs %d genome-wide rows" % (len(rows), len(genome)))

    for r in rows:
        key = (r["POP_GROUPING"], int(r["G"]))
        want = genome.get(key)
        got = [r["NUM_LOCI"], r["MEAN"], r["VAR"], r["STD_ERR"]]
        s.check("whole.%s.g%s" % key, want == got, "%s vs %s" % (want, got))


def case_recompute(s, adze, w):
    """Every window row equals a mean/variance computed here from the per-locus
    values the same run printed."""
    meta = gen_data.write_pair(os.path.join(w, "rec"), npops=3, nind=6,
                               nloci=24, nall=6, missing=0.05, seed=12,
                               chroms=3, spacing=1000)
    run(adze, w, ["--data", "rec.vcf", "--samples", "rec.samples",
                  "--max-g", "4", "--window-bp", "3000", "--step-bp", "3000",
                  "--full-richness", "--out-prefix", "rec"])

    names, per = read_fulldata(os.path.join(w, "rec.richness_fulldata"))
    coords = dict(zip(meta["names"], meta["coords"]))
    order = [n for n in names]

    bad = 0
    for r in read_windows(os.path.join(w, "rec.richness_windows")):
        lo, hi = int(r["START"]), int(r["END"])
        idx = [i for i, n in enumerate(order)
               if coords[n][0] == r["CHROM"] and lo <= coords[n][1] <= hi]
        vals = [per[(r["POP_GROUPING"], int(r["G"]))][i] for i in idx]

        if len(vals) != int(r["NUM_LOCI"]):
            bad += 1
            s.check("recompute.nloci", False,
                    "%s:%s-%s says %s loci, layout gives %d"
                    % (r["CHROM"], lo, hi, r["NUM_LOCI"], len(vals)))
            continue

        n = len(vals)
        mean = sum(vals) / n
        if not close(mean, float(r["MEAN"])):
            bad += 1
            s.check("recompute.mean", False,
                    "%s:%s-%s %s g%s: %.8g vs %s" % (r["CHROM"], lo, hi,
                    r["POP_GROUPING"], r["G"], mean, r["MEAN"]))
            continue
        if n > 1:
            var = sum((v - mean) ** 2 for v in vals) / (n - 1)
            se = (var / n) ** 0.5
            if not (close(var, float(r["VAR"]), VAR_RTOL) and
                    close(se, float(r["STD_ERR"]), VAR_RTOL)):
                bad += 1
                s.check("recompute.var", False,
                        "%s:%s-%s: var %.8g vs %s" % (r["CHROM"], lo, hi,
                        var, r["VAR"]))
    s.check("recompute.all", bad == 0, "%d rows disagreed" % bad)


def case_layout(s, adze, w):
    """Windows tile as asked, stay inside a chromosome, and overlap by half."""
    meta = gen_data.write_pair(os.path.join(w, "lay"), npops=3, nind=5,
                               nloci=18, nall=4, seed=13,
                               chroms=[10, 5, 3], spacing=1000)
    coords = dict(zip(meta["names"], meta["coords"]))

    run(adze, w, ["--data", "lay.vcf", "--samples", "lay.samples",
                  "--max-g", "3", "--window-bp", "4000", "--step-bp", "2000",
                  "--out-prefix", "lay"])
    rows = read_windows(os.path.join(w, "lay.richness_windows"))

    # One grouping's worth of rows at one g is the window list itself.
    wins = [(r["CHROM"], int(r["START"]), int(r["END"]), int(r["NUM_LOCI"]))
            for r in rows if r["POP_GROUPING"] == "POP0" and r["G"] == "2"]

    # 4 kb windows stepping 2 kb over 10, 5 and 3 loci spaced 1 kb apart.
    want = [("chr1", 1, 4000, 4), ("chr1", 2001, 6000, 4),
            ("chr1", 4001, 8000, 4), ("chr1", 6001, 10000, 4),
            ("chr1", 8001, 12000, 2),
            ("chr2", 1, 4000, 4), ("chr2", 2001, 6000, 3),
            ("chr2", 4001, 8000, 1),
            ("chr3", 1, 4000, 3), ("chr3", 2001, 6000, 1)]
    s.check("layout.windows", wins == want, "%s" % (wins,))

    # Every window lies inside one chromosome, and its locus count matches the
    # loci actually in its interval.
    for c, lo, hi, n in wins:
        inside = sum(1 for nm in meta["names"]
                     if coords[nm][0] == c and lo <= coords[nm][1] <= hi)
        s.check("layout.%s:%d" % (c, lo), inside == n,
                "%d loci in the interval, row says %d" % (inside, n))

    # With the step at half the width, an interior locus falls in two windows.
    counts = {}
    for c, lo, hi, _ in wins:
        for nm in meta["names"]:
            if coords[nm][0] == c and lo <= coords[nm][1] <= hi:
                counts[nm] = counts.get(nm, 0) + 1
    interior = [nm for nm in meta["names"]
                if coords[nm][0] == "chr1" and 2000 < coords[nm][1] < 9000]
    s.check("layout.overlap", all(counts[nm] == 2 for nm in interior),
            "%s" % {nm: counts[nm] for nm in interior})


def case_units(s, adze, w):
    """Evenly spaced loci: basepair and locus windows cover the same loci."""
    gen_data.write_pair(os.path.join(w, "un"), npops=3, nind=5, nloci=20,
                        nall=5, seed=14, chroms=2, spacing=1000)
    common = ["--data", "un.vcf", "--samples", "un.samples", "--max-g", "4"]
    run(adze, w, common + ["--window-bp", "5000", "--out-prefix", "un_bp"])
    run(adze, w, common + ["--window-loci", "5", "--out-prefix", "un_loci"])

    a = read_windows(os.path.join(w, "un_bp.richness_windows"))
    b = read_windows(os.path.join(w, "un_loci.richness_windows"))
    s.check("units.rowcount", len(a) == len(b),
            "%d vs %d rows" % (len(a), len(b)))

    # START and END differ by construction: a basepair window reports its own
    # bounds, a locus window the span of the loci it holds.
    keys = ["CHROM", "POP_GROUPING", "G", "NUM_LOCI", "MEAN", "VAR", "STD_ERR"]
    same = all([x[k] for k in keys] == [y[k] for k in keys]
               for x, y in zip(a, b))
    s.check("units.values", same, "basepair and locus windows disagree")


def case_minloci(s, adze, w):
    """--min-window-loci drops exactly the windows below the threshold."""
    gen_data.write_pair(os.path.join(w, "mn"), npops=3, nind=5, nloci=18,
                        nall=4, seed=15, chroms=[10, 5, 3], spacing=1000)
    common = ["--data", "mn.vcf", "--samples", "mn.samples", "--max-g", "3",
              "--window-bp", "4000", "--step-bp", "2000"]
    run(adze, w, common + ["--out-prefix", "mn_all"])
    run(adze, w, common + ["--min-window-loci", "3", "--out-prefix", "mn_3"])

    a = [r for r in read_windows(os.path.join(w, "mn_all.richness_windows"))
         if r["POP_GROUPING"] == "POP0" and r["G"] == "2"]
    b = [r for r in read_windows(os.path.join(w, "mn_3.richness_windows"))
         if r["POP_GROUPING"] == "POP0" and r["G"] == "2"]
    want = [r for r in a if int(r["NUM_LOCI"]) >= 3]
    s.check("minloci.kept",
            [(r["CHROM"], r["START"]) for r in b] ==
            [(r["CHROM"], r["START"]) for r in want],
            "kept %d, expected %d" % (len(b), len(want)))


def case_formats(s, adze, w):
    """VCF and STRUCTURE-plus-map window to the same output."""
    gen_data.write_pair(os.path.join(w, "fm"), npops=3, nind=5, nloci=18,
                        nall=5, missing=0.05, seed=16, chroms=2, spacing=1000)
    run(adze, w, ["--data", "fm.vcf", "--samples", "fm.samples",
                  "--max-g", "4", "--window-loci", "4", "--out-prefix", "fm_v"])
    run(adze, w, ["--data", "fm.stru", "--group-col", "2", "--loci-map",
                  "fm.map", "--max-g", "4", "--window-loci", "4",
                  "--out-prefix", "fm_s"])
    for stat in ("richness_windows", "private_windows"):
        a = open(os.path.join(w, "fm_v." + stat)).read()
        b = open(os.path.join(w, "fm_s." + stat)).read()
        s.check("formats.%s" % stat, a == b, "VCF and STRUCTURE output differ")


def case_threads(s, adze, w):
    """Window output is the same at any thread count."""
    gen_data.write_pair(os.path.join(w, "th"), npops=4, nind=6, nloci=40,
                        nall=6, seed=17, chroms=2, spacing=1000)
    base = ["--data", "th.vcf", "--samples", "th.samples", "--max-g", "6",
            "--window-loci", "5"]
    outs = []
    for n in ("1", "4"):
        r = run(adze, w, base + ["--threads", n, "--out-prefix", "th_" + n])
        if "no OpenMP support" in r.stderr and n != "1":
            print("SKIP  threads (build has no OpenMP)")
            return
        outs.append(open(os.path.join(w, "th_%s.richness_windows" % n)).read())
    s.check("threads.identical", outs[0] == outs[1],
            "1 and 4 threads differ")


def case_refusals(s, adze, w):
    """Inputs a window cannot be defined over are refused, not guessed at."""
    meta = gen_data.write_pair(os.path.join(w, "rf"), npops=3, nind=5,
                               nloci=12, nall=4, seed=18, chroms=2,
                               spacing=1000)

    lines = open(os.path.join(w, "rf.vcf")).read().splitlines()
    head = [l for l in lines if l.startswith("#")]
    body = [l for l in lines if not l.startswith("#")]

    swapped = body[:]
    swapped[2], swapped[4] = swapped[4], swapped[2]
    open(os.path.join(w, "rf_unsorted.vcf"), "w").write(
        "\n".join(head + swapped) + "\n")

    mixed = body[:2] + [body[-1]] + body[2:-1]
    open(os.path.join(w, "rf_mixed.vcf"), "w").write(
        "\n".join(head + mixed) + "\n")

    full = open(os.path.join(w, "rf.map")).read().splitlines()
    open(os.path.join(w, "rf_short.map"), "w").write(
        "\n".join(full[:5]) + "\n")

    for label, args, want in [
            ("unsorted", ["--data", "rf_unsorted.vcf", "--samples",
                          "rf.samples", "--window-bp", "3000"], 4),
            ("chrom-blocks", ["--data", "rf_mixed.vcf", "--samples",
                              "rf.samples", "--window-bp", "3000"], 4),
            ("no-coordinate", ["--data", "rf.stru", "--group-col", "2",
                               "--loci-map", "rf_short.map",
                               "--window-loci", "3"], 4),
            ("no-map", ["--data", "rf.stru", "--group-col", "2",
                        "--window-loci", "3"], 2)]:
        r = run(adze, w, args + ["--out-prefix", "rf_" + label])
        s.check("refusals.%s" % label, r.returncode == want,
                "exit %d, expected %d" % (r.returncode, want))
        s.check("refusals.%s.message" % label, "ERROR" in r.stderr,
                "no ERROR line on stderr")

    # An unsorted file is still fine when no window is asked for.
    r = run(adze, w, ["--data", "rf_unsorted.vcf", "--samples", "rf.samples",
                      "--out-prefix", "rf_ok"])
    s.check("refusals.unsorted-ok-without-windows", r.returncode == 0,
            "exit %d" % r.returncode)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate", required=True, help="adze binary to test")
    ap.add_argument("--workdir", default=os.path.join(HERE, "work-windows"))
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    shutil.rmtree(args.workdir, ignore_errors=True)
    os.makedirs(args.workdir)

    s = Suite(args.verbose)
    for case in (case_whole, case_recompute, case_layout, case_units,
                 case_minloci, case_formats, case_threads, case_refusals):
        case(s, args.candidate, args.workdir)

    if not s.failures and not args.keep:
        shutil.rmtree(args.workdir, ignore_errors=True)

    print("\n%d passed, %d failed" % (s.npass, len(s.failures)))
    if s.failures:
        print("scratch kept: %s" % args.workdir)
    return 1 if s.failures else 0


if __name__ == "__main__":
    sys.exit(main())
