#!/usr/bin/env python3
"""Hold a candidate build to the behaviour of another build of ADZE 2.x.

Where regress.py holds the numbers to version 1.0 -- over the rows 1.0 can
produce, with a tolerance on the summary columns -- this suite holds a
candidate to a *2.x reference*: same inputs, same flags, byte-identical
output. It exists for the streaming rewrite, where the arithmetic is meant to
be untouched and every difference is therefore a defect until it is one of
the decisions recorded in the plan.

    ./test/equivalence.py --candidate src/adze --reference /path/to/adze-pre

Both binaries run in their own scratch directory, and every file either wrote
is compared. Lines carrying a wall-clock duration are masked, since two runs
of the same binary differ there.

A decision that deliberately changes a file can be declared, so that it shows
up as a named exception rather than being quietly tolerated. The declaration
says exactly what the new file must be, so it is still an exact comparison:

    --fulldata-locus-major      the reference's _fulldata, with its label and
                                LOCUS columns swapped and its rows re-sorted
                                by locus (groupings keeping their order within
                                a locus), must equal the candidate's byte for
                                byte. See doc/streaming.md, decision 1.

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

TIME_LINE = re.compile(r"\(d:h:m:s[^)]*\)|\bseconds\b|\bs\)\s*$")

# Fixtures: shapes that have caught something before, plus a windowed one.
FIXTURES = [
    ("basic",     dict(npops=4, nind=8,  nloci=30,  nall=6,  missing=0.00, seed=1)),
    ("missing",   dict(npops=5, nind=10, nloci=40,  nall=8,  missing=0.25, seed=2)),
    ("biallelic", dict(npops=6, nind=12, nloci=120, nall=2,  missing=0.05, seed=3)),
    ("rich",      dict(npops=3, nind=15, nloci=25,  nall=40, missing=0.10, seed=4)),
    ("uneven",    dict(npops=4, nind=[4, 9, 15, 25], nloci=30, nall=7, missing=0.15, seed=7)),
    ("haploid",   dict(npops=3, nind=12, nloci=20,  nall=6,  missing=0.05, seed=8, ploidy=1)),
    ("twochrom",  dict(npops=4, nind=10, nloci=60,  nall=5,  missing=0.05, seed=10)),
]

# Flag sets, applied to whichever encodings they make sense for.
STRU = ["--data", "f.stru", "--group-col", "2"]
VCF = ["--data", "f.vcf", "--samples", "f.samples"]
VCFGZ = ["--data", "f.vcf.gz", "--samples", "f.samples"]

CASES = [
    ("stru.default",    STRU, []),
    ("stru.notol",      STRU, ["--tolerance", "1"]),
    ("stru.ladder8",    STRU, ["--tolerance", "1", "--max-g", "8"]),
    ("stru.atg",        STRU, ["--tolerance", "1", "--at-g", "5"]),
    ("stru.legacy",     STRU, ["--tolerance", "1", "--legacy"]),
    ("stru.full",       STRU, ["--tolerance", "1", "--full-richness", "--full-private"]),
    ("stru.tuples",     STRU, ["--tolerance", "1", "--combinations", "--tuples-k", "1,2",
                               "--full-tuples"]),
    ("stru.richonly",   STRU, ["--tolerance", "1", "--stat", "richness"]),
    ("stru.pops",       STRU, ["--tolerance", "1", "--pops", "POP0,POP1"]),
    ("stru.windows",    STRU, ["--tolerance", "1", "--loci-map", "f.map",
                               "--window-loci", "5", "--step-loci", "5",
                               "--min-window-loci", "2", "--full-richness"]),
    ("vcf.default",     VCF,  []),
    ("vcf.notol",       VCF,  ["--tolerance", "1"]),
    ("vcf.atg",         VCF,  ["--tolerance", "1", "--at-g", "4", "--full-private"]),
    ("vcf.tuples",      VCF,  ["--tolerance", "1", "--combinations", "--tuples-k", "2"]),
    ("vcf.windowsbp",   VCF,  ["--tolerance", "1", "--window-bp", "4000", "--step-bp", "2000",
                               "--min-window-loci", "2"]),
    ("vcfgz.default",   VCFGZ, ["--tolerance", "1"]),
    ("vcf.dryrun",      VCF,  ["--dry-run"]),
    # The dry run is served from the scan rather than from a loaded dataset,
    # so each thing it reports needs a case: dimensions, the window layout,
    # the tuple count, the clamped --at-g, and a filtered set of groupings.
    ("stru.dryrun",     STRU, ["--dry-run"]),
    ("stru.dryrun.pops", STRU, ["--pops", "POP0,POP1", "--dry-run"]),
    ("stru.dryrun.win", STRU, ["--tolerance", "1", "--loci-map", "f.map",
                               "--window-loci", "5", "--step-loci", "5", "--dry-run"]),
    ("vcf.dryrun.win",  VCF,  ["--window-bp", "4000", "--step-bp", "2000",
                               "--min-window-loci", "2", "--dry-run"]),
    ("vcf.dryrun.tuples", VCF, ["--combinations", "--tuples-k", "1,2", "--dry-run"]),
    ("vcf.dryrun.atg",  VCF,  ["--at-g", "max", "--dry-run"]),
    ("vcf.dryrun.notol", VCF, ["--tolerance", "1", "--dry-run"]),
]


def fixture(workdir, name, kw):
    d = os.path.join(workdir, "data", name)
    os.makedirs(d, exist_ok=True)
    stamp = os.path.join(d, "f.stru")
    if not os.path.exists(stamp):
        chroms = 2 if name == "twochrom" else 1
        gen_data.write_pair(os.path.join(d, "f"), chroms=chroms, **kw)
    return d


def masked(path):
    with open(path, "rb") as fh:
        raw = fh.read()
    try:
        text = raw.decode()
    except UnicodeDecodeError:
        return raw
    lines = [l for l in text.splitlines() if not TIME_LINE.search(l)]
    if path.endswith("_summary") or path.endswith(".summary.txt"):
        lines = [l for l in lines if l.strip()]
    return "\n".join(lines).encode()


def first_difference(a, b):
    ra, rb = a.decode(errors="replace").splitlines(), b.decode(errors="replace").splitlines()
    for i, (x, y) in enumerate(zip(ra, rb)):
        if x != y:
            return "line %d: ref %r vs cnd %r" % (i + 1, x[:64], y[:64])
    if len(ra) != len(rb):
        longer, which = (ra, "reference") if len(ra) > len(rb) else (rb, "candidate")
        return "%s has %d more lines, first %r" % (
            which, abs(len(ra) - len(rb)), longer[min(len(ra), len(rb))][:64])
    return "differ in trailing bytes only"


def outputs(directory):
    skip = {"f.stru", "f.vcf", "f.vcf.gz", "f.samples", "f.map",
            "console.out", "console.err"}
    return sorted(f for f in os.listdir(directory) if f not in skip)


def as_locus_major(text):
    """The reference's _fulldata, rewritten as the locus-major layout.

    Decision 1 of the streaming rewrite: rows follow the sweep, so they come
    out ordered by locus with a locus's groupings together, and LOCUS leads
    the row because it is the sort key. Nothing else about the file changes,
    so the expected content is computable from the reference's own and the
    comparison stays exact — a changed value still fails.
    """
    lines = text.splitlines()
    if not lines:
        return text
    head = lines[0].split("\t")
    if len(head) < 3 or head[1] != "LOCUS":
        return text                      # not a file this decision covers

    rows = [l.split("\t") for l in lines[1:] if l.strip()]
    loci, groups = [], []
    for r in rows:                       # first-appearance order, both axes
        if r[1] not in loci:
            loci.append(r[1])
        if r[0] not in groups:
            groups.append(r[0])
    at = {(r[0], r[1]): r for r in rows}

    out = ["\t".join([head[1], head[0]] + head[2:])]
    for locus in loci:
        for g in groups:
            r = at.get((g, locus))
            if r is not None:
                out.append("\t".join([r[1], r[0]] + r[2:]))
    return "\n".join(out) + "\n"


def compare(ref_dir, cnd_dir, reordered_fulldata):
    ref_files, cnd_files = outputs(ref_dir), outputs(cnd_dir)
    if ref_files != cnd_files:
        only_ref = [f for f in ref_files if f not in cnd_files]
        only_cnd = [f for f in cnd_files if f not in ref_files]
        return "file set differs (reference only: %s; candidate only: %s)" % (
            only_ref or "-", only_cnd or "-")

    for f in ref_files:
        a, b = masked(os.path.join(ref_dir, f)), masked(os.path.join(cnd_dir, f))
        if a == b:
            continue
        if reordered_fulldata and f.endswith("_fulldata"):
            # A declared decision: the expected file is computed from the
            # reference's own content, so this is still an exact comparison.
            expected = as_locus_major(a.decode()).rstrip("\n").encode()
            if expected == b.rstrip(b"\n"):
                continue
            return "%s (not the declared locus-major rewrite: %s)" % (
                f, first_difference(expected, b.rstrip(b"\n")))
        return "%s (%s)" % (f, first_difference(a, b))

    return None


def run(binary, workdir, argv):
    proc = subprocess.run([os.path.abspath(binary)] + argv, cwd=workdir,
                          capture_output=True, text=True)
    with open(os.path.join(workdir, "console.out"), "w", newline="\n") as fh:
        fh.write(proc.stdout)
    with open(os.path.join(workdir, "console.err"), "w", newline="\n") as fh:
        fh.write(proc.stderr)
    return proc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--workdir", default=os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "work-equivalence"))
    ap.add_argument("--fulldata-locus-major", action="store_true",
                    help="expect _fulldata in the locus-major layout (doc/streaming.md, decision 1)")
    ap.add_argument("--only", default=None, help="substring filter on case names")
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    shutil.rmtree(args.workdir, ignore_errors=True)
    os.makedirs(args.workdir, exist_ok=True)

    npass = nfail = 0
    failures = []
    for fname, kw in FIXTURES:
        data = fixture(args.workdir, fname, kw)
        for case, base, extra in CASES:
            name = "%s.%s" % (fname, case)
            if args.only and args.only not in name:
                continue
            dirs = {}
            for side, binary in (("ref", args.reference), ("cnd", args.candidate)):
                d = os.path.join(args.workdir, name, side)
                os.makedirs(d, exist_ok=True)
                for f in os.listdir(data):
                    shutil.copy(os.path.join(data, f), d)
                dirs[side] = d
            argv = base + extra + ["--out-prefix", "o"]
            r = run(args.reference, dirs["ref"], argv)
            c = run(args.candidate, dirs["cnd"], argv)

            if r.returncode != c.returncode:
                why = "exit status %d vs %d" % (r.returncode, c.returncode)
            elif r.stdout != c.stdout and masked_text(r.stdout) != masked_text(c.stdout):
                why = "stdout differs: %s" % first_difference(
                    masked_text(r.stdout).encode(), masked_text(c.stdout).encode())
            else:
                why = compare(dirs["ref"], dirs["cnd"], args.fulldata_locus_major)

            if why:
                nfail += 1
                failures.append((name, why, dirs["cnd"]))
                print("FAIL  %-34s %s" % (name, why))
            else:
                npass += 1
                if args.verbose:
                    print("PASS  %-34s %d files" % (name, len(outputs(dirs["ref"]))))

    if not failures and not args.keep:
        shutil.rmtree(args.workdir, ignore_errors=True)

    print("\n%d passed, %d failed" % (npass, nfail))
    if failures:
        print("scratch kept: %s" % args.workdir)
    return 1 if failures else 0


def masked_text(text):
    return "\n".join(l for l in text.splitlines() if not TIME_LINE.search(l))


if __name__ == "__main__":
    sys.exit(main())
