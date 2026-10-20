#!/usr/bin/env python3
"""Haplotype blocks, against an implementation that shares no code with ADZE.

The fixtures carry their own haplotype matrix (gen_data.write_phased writes
<prefix>.haps before it writes a single record), so everything here is checked
against haplotypes assembled from that matrix rather than against ADZE's own
account of them.

What each group of cases holds:

  layout     block boundaries, spans and SNP counts, for --hap-snps and
             --hap-bp, on one and two chromosomes
  callable   the gene copies each block can score, under the rule that one
             uncalled genotype makes that copy's whole haplotype uncalled
  filter     which blocks --tolerance drops, and the MAX_G that follows

The scan's own numbers come from ADZE_DUMP_SCAN, which writes one row per
block and grouping; nothing here reads a statistic, so these cases isolate
the counting pass from the arithmetic that consumes it.

    ./test/haplotypes.py --candidate src/adze [--workdir DIR] [-v]

Exit status is 0 only if every case passes.
"""
import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_data


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


def run(binary, workdir, args, env=None):
    return subprocess.run([os.path.abspath(binary)] + args, cwd=workdir,
                          capture_output=True, text=True,
                          env=dict(os.environ, **env) if env else None)


def scan_rows(path):
    """{(block index, grouping): (callable copies, missing copies)}."""
    out = {}
    with open(path) as fh:
        head = fh.readline()
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            out[(int(f[0]), f[1])] = (int(f[2]), int(f[3]))
    return out


def records(vcf):
    """[(chromosome, position)] in file order, and the sample names."""
    coords, samples = [], []
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]
                continue
            coords.append((f[0], int(f[1])))
    return coords, samples


def blocks_by_snps(coords, n, minsnps=1):
    """Block layout as the program should produce it: runs of n SNPs."""
    out, cur = [], []
    for i, (chrom, pos) in enumerate(coords):
        if cur and (chrom != coords[cur[0]][0] or len(cur) >= n):
            out.append(cur)
            cur = []
        cur.append(i)
    if cur:
        out.append(cur)
    return [b for b in out if len(b) >= minsnps]


def blocks_by_bp(coords, width, minsnps=1):
    """The same, as fixed spans anchored at position 1."""
    out, cur, key = [], [], None
    for i, (chrom, pos) in enumerate(coords):
        k = (chrom, (pos - 1) // width)
        if cur and k != key:
            out.append(cur)
            cur = []
        cur.append(i)
        key = k
    if cur:
        out.append(cur)
    return [b for b in out if len(b) >= minsnps]


def callable_copies(haps, groups, block):
    """Copies whose haplotype is callable across `block`, per grouping.

    One uncalled genotype anywhere in the block makes that copy uncalled for
    the whole block, so this is an all-or-nothing test over the block's SNPs.
    """
    out = {}
    for (sample, copy), seq in haps.items():
        if all(seq[i] is not None for i in block):
            g = groups[sample]
            out[g] = out.get(g, 0) + 1
    return out


def groups_of(samples_file):
    out = {}
    for line in open(samples_file):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) >= 2:
            out[f[0]] = f[1]
    return out


def case_layout(s, adze, w):
    """Block boundaries, spans and SNP counts come out as laid down."""
    for name, kw in gen_data.PHASED:
        if "unphased" in kw:
            continue
        meta = gen_data.write_phased(os.path.join(w, name), **kw)
        coords, _ = records(os.path.join(w, meta["file"]))

        for flag, width, expected in (("--hap-snps", 5, blocks_by_snps(coords, 5)),
                                      ("--hap-bp", 5000, blocks_by_bp(coords, 5000))):
            dump = os.path.join(w, "%s%s.scan" % (name, flag[2:]))
            r = run(adze, w, ["--data", meta["file"], "--samples", meta["samples"],
                              "--tolerance", "1", flag, str(width), "--dry-run"],
                    env={"ADZE_DUMP_SCAN": dump})
            s.check("layout.%s.%s.runs" % (name, flag[6:]), r.returncode == 0,
                    "exit %d: %s" % (r.returncode, r.stderr.strip()[:120]))
            if r.returncode != 0:
                continue

            said = "collapsed into %d haplotype block" % len(expected)
            s.check("layout.%s.%s.count" % (name, flag[6:]),
                    said in r.stdout + r.stderr,
                    "expected %d blocks; run said: %s" % (
                        len(expected),
                        [l for l in (r.stdout + r.stderr).splitlines()
                         if "collapsed" in l][:1]))

            rows = scan_rows(dump)
            s.check("layout.%s.%s.rows" % (name, flag[6:]),
                    len({b for b, _ in rows}) == len(expected),
                    "scan wrote %d blocks, layout gives %d" % (
                        len({b for b, _ in rows}), len(expected)))


def case_callable(s, adze, w):
    """Every block's callable copies, against the fixture's own haplotypes."""
    for name, kw in gen_data.PHASED:
        if "unphased" in kw:
            continue
        meta = gen_data.write_phased(os.path.join(w, name), **kw)
        coords, _ = records(os.path.join(w, meta["file"]))
        _, haps = gen_data.read_haps(os.path.join(w, meta["haps"]))
        groups = groups_of(os.path.join(w, meta["samples"]))

        for width in (2, 5):
            dump = os.path.join(w, "%s_c%d.scan" % (name, width))
            r = run(adze, w, ["--data", meta["file"], "--samples", meta["samples"],
                              "--tolerance", "1", "--hap-snps", str(width),
                              "--dry-run"],
                    env={"ADZE_DUMP_SCAN": dump})
            if r.returncode != 0:
                s.check("callable.%s.%d.runs" % (name, width), False,
                        "exit %d" % r.returncode)
                continue

            rows = scan_rows(dump)
            want = {}
            for b, block in enumerate(blocks_by_snps(coords, width)):
                for g, n in callable_copies(haps, groups, block).items():
                    want[(b, g)] = n

            bad = [(k, want.get(k, 0), rows.get(k, (None,))[0])
                   for k in set(list(want) + list(rows))
                   if want.get(k, 0) != rows.get(k, (0,))[0]]
            s.check("callable.%s.snps%d" % (name, width), not bad,
                    "%d of %d block-grouping cells differ, e.g. %s" % (
                        len(bad), max(len(want), len(rows)), bad[:3]))


def case_filter(s, adze, w):
    """Which blocks --tolerance drops, and the ceiling that follows."""
    name, kw = "hap_missing", dict(npops=4, nind=6, nloci=24, nall=2,
                                   missing=0.10, seed=12)
    meta = gen_data.write_phased(os.path.join(w, name), **kw)
    coords, _ = records(os.path.join(w, meta["file"]))
    _, haps = gen_data.read_haps(os.path.join(w, meta["haps"]))
    groups = groups_of(os.path.join(w, meta["samples"]))

    copies = {}
    for (sample, _), _seq in haps.items():
        g = groups[sample]
        copies[g] = copies.get(g, 0) + 1

    for tol in ("1", "0.5", "0.2"):
        r = run(adze, w, ["--data", meta["file"], "--samples", meta["samples"],
                          "--hap-snps", "4", "--tolerance", tol,
                          "--out-prefix", "f%s" % tol.replace(".", "")])
        if r.returncode not in (0, 4):
            s.check("filter.tol%s.runs" % tol, False, "exit %d" % r.returncode)
            continue

        blocks = blocks_by_snps(coords, 4)
        kept = []
        for block in blocks:
            have = callable_copies(haps, groups, block)
            worst = max((copies[g] - have.get(g, 0)) / float(copies[g])
                        for g in copies)
            if worst <= float(tol):
                kept.append(block)

        if not kept:
            s.check("filter.tol%s.refused" % tol, r.returncode == 4,
                    "no block should survive, but the run exited %d" % r.returncode)
            continue

        ceiling = min(min(callable_copies(haps, groups, b).get(g, 0)
                          for b in kept) for g in copies)
        said = "using %d, the largest" % ceiling
        s.check("filter.tol%s.maxg" % tol,
                said in r.stdout + r.stderr or ceiling == 0,
                "expected MAX_G %d over %d surviving blocks; run said: %s" % (
                    ceiling, len(kept),
                    [l for l in (r.stdout + r.stderr).splitlines()
                     if "MAX_G" in l][:1]))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--workdir", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "work-haplotypes"))
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    s = Suite(args.verbose)
    for case in (case_layout, case_callable, case_filter):
        case(s, args.candidate, args.workdir)

    print("\n%d passed, %d failed" % (s.passed, len(s.failed)))
    return 1 if s.failed else 0


if __name__ == "__main__":
    sys.exit(main())
