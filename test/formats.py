#!/usr/bin/env python3
"""Check that the STRUCTURE, VCF and compressed-VCF readers agree.

The same genotypes are written in all three encodings (test/gen_data.py) and
run through the same binary. Two comparisons are made:

  vcf vs vcf.gz       byte-identical: same reader, same input, only the
                      container differs
  structure vs vcf    byte-identical: both readers key allele slots the same
                      way (grouping, then position within it, then gene copy),
                      so Nji comes out in the same layout and every statistic
                      sums its terms in the same order
  lf vs crlf          byte-identical: a file with Windows line endings is the
                      same data, so the trailing carriage return must not
                      reach an allele label

Byte equality is the point. Two readers that merely agreed to within rounding
would leave it open whether a later change had altered the numbers or only the
summation order; with the ordering rule shared, any difference at all is a
real one.

    ./test/formats.py --candidate src/adze

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

CASES = [
    ("basic",     dict(npops=4, nind=8,  nloci=30, nall=6, missing=0.00, seed=1)),
    ("missing",   dict(npops=5, nind=10, nloci=25, nall=8, missing=0.25, seed=2)),
    ("biallelic", dict(npops=6, nind=12, nloci=60, nall=2, missing=0.05, seed=3)),
    ("uneven",    dict(npops=4, nind=[4, 9, 15, 25], nloci=20, nall=7,
                       missing=0.15, seed=7)),
    ("haploid",   dict(npops=3, nind=12, nloci=20, nall=6, missing=0.05,
                       seed=8, ploidy=1)),
]


def run(binary, workdir, args):
    r = subprocess.run([os.path.abspath(binary)] + args, cwd=workdir,
                       capture_output=True, text=True)
    return r


def same_bytes(a, b):
    """None if the two files are identical, else a description of the first
    line that differs."""
    la = open(a).read().splitlines()
    lb = open(b).read().splitlines()
    if len(la) != len(lb):
        return "%d lines vs %d lines" % (len(la), len(lb))
    for n, (x, y) in enumerate(zip(la, lb), start=1):
        if x != y:
            return "line %d: %r vs %r" % (n, x[:60], y[:60])
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate", required=True, help="adze binary to test")
    ap.add_argument("--workdir", default=os.path.join(HERE, "work-formats"))
    ap.add_argument("--keep", action="store_true", help="keep scratch files")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    shutil.rmtree(args.workdir, ignore_errors=True)
    os.makedirs(args.workdir)

    npass = nfail = 0
    failures = []

    for name, kw in CASES:
        meta = gen_data.write_pair(os.path.join(args.workdir, name), **kw)
        # --tolerance 1 keeps every locus: this suite compares the readers
        # against each other, so the locus set must not depend on the default.
        common = ["--max-g", "6", "--full-richness", "--full-private",
                  "--tolerance", "1"]

        outs = {}
        for tag, source in [("stru", ["--data", "%s.stru" % name,
                                      "--group-col", "2"]),
                            ("vcf", ["--data", "%s.vcf" % name,
                                     "--samples", "%s.samples" % name]),
                            ("gz", ["--data", "%s.vcf.gz" % name,
                                    "--samples", "%s.samples" % name])]:
            r = run(args.candidate, args.workdir,
                    source + common + ["--out-prefix", "%s_%s" % (name, tag)])
            outs[tag] = r
            if r.returncode != 0 and tag != "gz":
                failures.append(("%s.%s" % (name, tag),
                                 "exit %d: %s" % (r.returncode,
                                                  r.stderr.strip().splitlines()[-1:])))

        def path(tag, stat):
            return os.path.join(args.workdir, "%s_%s.%s" % (name, tag, stat))

        # vcf vs vcf.gz: identical bytes, unless this build has no zlib
        if outs["gz"].returncode == 0:
            for stat in ("richness", "private", "richness_fulldata"):
                same = (open(path("vcf", stat), "rb").read() ==
                        open(path("gz", stat), "rb").read())
                label = "%s.vcf==gz.%s" % (name, stat)
                if same:
                    npass += 1
                    if args.verbose:
                        print("PASS  %s" % label)
                else:
                    nfail += 1
                    failures.append((label, "compressed and plain differ"))
        elif "no zlib support" in outs["gz"].stderr:
            print("SKIP  %s.vcf.gz (build has no zlib)" % name)
        else:
            nfail += 1
            failures.append(("%s.gz" % name, outs["gz"].stderr.strip()[-200:]))

        # CRLF input: the same genotypes with Windows line endings. Without
        # the reader stripping the carriage return, the last allele on each
        # row interns as its own label and every count at that locus shifts --
        # so this is checked on every platform, not only where CRLF is native.
        crlf_in = os.path.join(args.workdir, "%s_crlf.stru" % name)
        with open(os.path.join(args.workdir, "%s.stru" % name), "rb") as fh:
            raw = fh.read()
        with open(crlf_in, "wb") as fh:
            fh.write(raw.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
        outs["crlf"] = run(args.candidate, args.workdir,
                           ["--data", "%s_crlf.stru" % name, "--group-col", "2"]
                           + common + ["--out-prefix", "%s_crlf" % name])
        for stat in ("richness", "private", "richness_fulldata"):
            label = "%s.lf==crlf.%s" % (name, stat)
            why = same_bytes(path("stru", stat), path("crlf", stat))
            if outs["crlf"].returncode != 0:
                why = "exit %d: %s" % (outs["crlf"].returncode,
                                       outs["crlf"].stderr.strip()[-160:])
            if why is None:
                npass += 1
                if args.verbose:
                    print("PASS  %s" % label)
            else:
                nfail += 1
                failures.append((label, why))

        # structure vs vcf: identical output for identical genotypes
        for stat in ("richness", "private", "richness_fulldata",
                     "private_fulldata"):
            label = "%s.stru==vcf.%s" % (name, stat)
            why = same_bytes(path("stru", stat), path("vcf", stat))
            if why is None:
                npass += 1
                if args.verbose:
                    print("PASS  %s" % label)
            else:
                nfail += 1
                failures.append((label, why))

    if not failures and not args.keep:
        shutil.rmtree(args.workdir, ignore_errors=True)

    print("\n%d passed, %d failed" % (npass, nfail))
    for label, why in failures:
        print("  FAIL %s: %s" % (label, why))
    if failures:
        print("\nscratch kept: %s" % args.workdir)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
