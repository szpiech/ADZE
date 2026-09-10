#!/usr/bin/env python3
"""Synthetic datasets for the ADZE regression, format and benchmark suites.

STRUCTURE format: one row per haplotype (a diploid individual occupies two
consecutive rows). Row 0 holds the locus names; each data row is
<ind_id> <pop_label> <allele>... with missing alleles written as -9.

write_pair() emits the same genotypes as both STRUCTURE and VCF, plus the
sample-to-grouping map VCF input needs, which is what test/formats.py uses to
check that the two readers agree.

Uses only the standard library so it runs anywhere the compiler does.
"""
import argparse
import gzip
import os
import random


def write_dataset(path, npops, nind, nloci, nall, missing=0.0, seed=1,
                  all_missing_loci=(), ploidy=2, allele_base=100):
    """Write one dataset and return the parameters ADZE needs to read it.

    nind may be an int (same sample size everywhere) or a per-population list,
    which is how unequal Nj across groups gets exercised.
    """
    rng = random.Random(seed)
    if isinstance(nind, int):
        nind = [nind] * npops
    assert len(nind) == npops
    with open(path, "w") as fh:
        fh.write(" ".join("L%d" % l for l in range(nloci)) + "\n")
        for p in range(npops):
            # Population-specific allele pool, overlapping but not identical, so
            # that private alleles and private alleles of k-tuples are non-zero.
            pool = [allele_base + a for a in range(nall)]
            rng.shuffle(pool)
            pool = pool[: max(2, int(round(0.7 * nall)))] + pool[:1] * 2
            for i in range(nind[p]):
                for _ in range(ploidy):
                    row = []
                    for l in range(nloci):
                        if l in all_missing_loci or rng.random() < missing:
                            row.append("-9")
                        else:
                            row.append(str(rng.choice(pool)))
                    fh.write("ind%d_%d POP%d %s\n" % (p, i, p, " ".join(row)))
    return {"dlines": ploidy * sum(nind), "loci": nloci,
            "nd_rows": 1, "nd_cols": 2, "sort_by": 2}


def _genotypes(npops, nind, nloci, nall, missing=0.0, seed=1,
                all_missing_loci=(), ploidy=2, allele_base=100):
    """Genotype matrix as [pop][individual][locus] -> tuple of allele codes.

    None is a missing gene copy. Allele codes are allele_base + index, so the
    index is recoverable for the VCF encoding.
    """
    rng = random.Random(seed)
    if isinstance(nind, int):
        nind = [nind] * npops
    assert len(nind) == npops
    out = []
    for p in range(npops):
        # Population-specific allele pool, overlapping but not identical, so
        # that private alleles and private alleles of k-tuples are non-zero.
        pool = [allele_base + a for a in range(nall)]
        rng.shuffle(pool)
        pool = pool[: max(2, int(round(0.7 * nall)))] + pool[:1] * 2
        inds = []
        for _ in range(nind[p]):
            copies = []
            for _ in range(ploidy):
                row = []
                for l in range(nloci):
                    if l in all_missing_loci or rng.random() < missing:
                        row.append(None)
                    else:
                        row.append(rng.choice(pool))
                copies.append(row)
            inds.append(copies)
        out.append(inds)
    return out, nind


def write_pair(prefix, nloci, allele_base=100, gzip_vcf=True, **kw):
    """Write the same genotypes as <prefix>.stru, <prefix>.vcf[.gz] and
    <prefix>.samples, and return the metadata for the STRUCTURE file.

    Locus names match across the two encodings (the VCF carries them in ID), so
    the *_fulldata columns line up and the outputs can be compared directly.
    """
    geno, nind = _genotypes(nloci=nloci, allele_base=allele_base, **kw)
    names = ["L%d" % l for l in range(nloci)]
    ploidy = len(geno[0][0])

    with open(prefix + ".stru", "w") as fh:
        fh.write(" ".join(names) + "\n")
        for p, inds in enumerate(geno):
            for i, copies in enumerate(inds):
                for row in copies:
                    fh.write("ind%d_%d POP%d %s\n" % (
                        p, i, p, " ".join("-9" if a is None else str(a) for a in row)))

    with open(prefix + ".samples", "w") as fh:
        fh.write("# sample\tgrouping\n")
        for p, inds in enumerate(geno):
            for i in range(len(inds)):
                fh.write("ind%d_%d\tPOP%d\n" % (p, i, p))

    samples = ["ind%d_%d" % (p, i) for p, inds in enumerate(geno)
               for i in range(len(inds))]
    nalt = max(a for inds in geno for c in inds for row in c
               for a in row if a is not None) - allele_base

    lines = ["##fileformat=VCFv4.2",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                              "FILTER", "INFO", "FORMAT"] + samples)]
    for l in range(nloci):
        # REF is allele index 0, so the ALT list covers indices 1..nalt.
        alt = ",".join("A" * (i + 1) for i in range(nalt)) if nalt else "."
        calls = []
        for p, inds in enumerate(geno):
            for copies in inds:
                calls.append("/".join(
                    "." if c[l] is None else str(c[l] - allele_base)
                    for c in copies))
        lines.append("\t".join(["chr1", str(l + 1), names[l], "T", alt, ".",
                                 "PASS", ".", "GT"] + calls))
    text = "\n".join(lines) + "\n"

    with open(prefix + ".vcf", "w") as fh:
        fh.write(text)
    if gzip_vcf:
        with gzip.open(prefix + ".vcf.gz", "wt") as fh:
            fh.write(text)

    return {"dlines": ploidy * sum(nind), "loci": nloci,
            "nd_rows": 1, "nd_cols": 2, "sort_by": 2,
            "file": os.path.basename(prefix) + ".stru",
            "samples": os.path.basename(prefix) + ".samples"}


def write_paramfile(path, meta, dfile, prefix, g=6, comb=0, k="", tol=1,
                    pp=0, full_r=0, full_p=0, full_c=0, missing="-9"):
    lines = ["MAX_G %d" % g,
             "DATA_LINES %d" % meta["dlines"],
             "LOCI %d" % meta["loci"],
             "NON_DATA_ROWS %d" % meta["nd_rows"],
             "NON_DATA_COLS %d" % meta["nd_cols"],
             "GROUP_BY_COL %d" % meta["sort_by"],
             "DATA_FILE %s" % dfile,
             "R_OUT %s_r" % prefix,
             "P_OUT %s_p" % prefix,
             "COMB %d" % comb,
             "C_OUT %s_c" % prefix,
             "MISSING %s" % missing,
             "TOLERANCE %s" % tol,
             "FULL_R %d" % full_r,
             "FULL_P %d" % full_p,
             "FULL_C %d" % full_c,
             "PRINT_PROGRESS %d" % pp]
    if comb:
        lines.append("K_RANGE %s" % k)
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


# The dataset matrix shared by the regression suite. Each entry is
# (name, kwargs for write_dataset) and deliberately covers the shapes where the
# 1.0 code paths differ: allele counts, missingness, one locus, one group,
# unequal sample sizes, and a locus with no observed alleles at all.
DATASETS = [
    ("basic",        dict(npops=4, nind=8,  nloci=30,  nall=6,  missing=0.00, seed=1)),
    ("missing",      dict(npops=5, nind=10, nloci=40,  nall=8,  missing=0.25, seed=2)),
    ("biallelic",    dict(npops=6, nind=12, nloci=200, nall=2,  missing=0.05, seed=3)),
    ("many_alleles", dict(npops=3, nind=15, nloci=25,  nall=40, missing=0.10, seed=4)),
    ("one_locus",    dict(npops=3, nind=10, nloci=1,   nall=5,  missing=0.00, seed=5)),
    ("one_group",    dict(npops=1, nind=20, nloci=20,  nall=6,  missing=0.10, seed=6)),
    ("uneven",       dict(npops=4, nind=[4, 9, 15, 25], nloci=30, nall=7, missing=0.15, seed=7)),
    ("haploid",      dict(npops=3, nind=12, nloci=20, nall=6, missing=0.05, seed=8, ploidy=1)),
    ("all_missing_locus",
     dict(npops=4, nind=10, nloci=6, nall=5, missing=0.05, seed=9, all_missing_loci=(3,))),
]


def build_all(outdir):
    os.makedirs(outdir, exist_ok=True)
    made = {}
    for name, kw in DATASETS:
        path = os.path.join(outdir, "%s.stru" % name)
        made[name] = write_dataset(path, **kw)
        made[name]["file"] = os.path.basename(path)
    return made


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("outdir")
    args = ap.parse_args()
    for name, meta in sorted(build_all(args.outdir).items()):
        print("%-20s %s" % (name, meta))
