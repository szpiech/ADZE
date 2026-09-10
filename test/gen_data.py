#!/usr/bin/env python3
"""Synthetic STRUCTURE-format datasets for the ADZE regression and benchmark suites.

One row per haplotype (a diploid individual occupies two consecutive rows).
Row 0 holds the locus names; each data row is  <ind_id> <pop_label> <allele>...
Missing alleles are written as -9.

Uses only the standard library so it runs anywhere the compiler does.
"""
import argparse
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
