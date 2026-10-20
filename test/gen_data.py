#!/usr/bin/env python3
"""Synthetic datasets for the ADZE regression, format and benchmark suites.

STRUCTURE format: one row per haplotype (a diploid individual occupies two
consecutive rows). Row 0 holds the locus names; each data row is
<ind_id> <pop_label> <allele>... with missing alleles written as -9.

write_pair() emits the same genotypes as both STRUCTURE and VCF, plus the
sample-to-grouping map VCF input needs, which is what test/formats.py uses to
check that the two readers agree.

Uses only the standard library so it runs anywhere the compiler does.

Every file is written with newline="\n", so a fixture is the same bytes on
every platform; Python's text mode would otherwise write CRLF on Windows. Both
builds read either ending correctly -- that is checked in test/formats.py --
but a fixture that changes with the platform makes a byte comparison harder to
reason about than it needs to be.
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
    with open(path, "w", newline="\n") as fh:
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


def place_loci(nloci, chroms=1, spacing=1000):
    """Coordinates for nloci loci, dealt out over `chroms` chromosomes.

    chroms may be an int (split as evenly as possible) or a list of per-
    chromosome locus counts. Positions start at `spacing` and step by it, so a
    window of k*spacing basepairs holds exactly k loci and basepair and
    locus-count windows can be predicted by hand.
    """
    if isinstance(chroms, int):
        base, extra = divmod(nloci, chroms)
        counts = [base + (1 if c < extra else 0) for c in range(chroms)]
    else:
        counts = list(chroms)
        assert sum(counts) == nloci, "locus counts must sum to nloci"
    out = []
    for c, n in enumerate(counts):
        for i in range(n):
            out.append(("chr%d" % (c + 1), spacing * (i + 1)))
    return out


def write_pair(prefix, nloci, allele_base=100, gzip_vcf=True, chroms=1,
               spacing=1000, **kw):
    """Write the same genotypes as <prefix>.stru, <prefix>.vcf[.gz],
    <prefix>.samples and <prefix>.map, and return the STRUCTURE metadata.

    Locus names match across the two encodings (the VCF carries them in ID), so
    the *_fulldata columns line up and the outputs can be compared directly.
    The coordinates in the VCF and in the locus map are the same, so a windowed
    run gives the same answer from either encoding.
    """
    geno, nind = _genotypes(nloci=nloci, allele_base=allele_base, **kw)
    names = ["L%d" % l for l in range(nloci)]
    ploidy = len(geno[0][0])
    coords = place_loci(nloci, chroms, spacing)

    with open(prefix + ".stru", "w", newline="\n") as fh:
        fh.write(" ".join(names) + "\n")
        for p, inds in enumerate(geno):
            for i, copies in enumerate(inds):
                for row in copies:
                    fh.write("ind%d_%d POP%d %s\n" % (
                        p, i, p, " ".join("-9" if a is None else str(a) for a in row)))

    with open(prefix + ".samples", "w", newline="\n") as fh:
        fh.write("# sample\tgrouping\n")
        for p, inds in enumerate(geno):
            for i in range(len(inds)):
                fh.write("ind%d_%d\tPOP%d\n" % (p, i, p))

    with open(prefix + ".map", "w", newline="\n") as fh:
        fh.write("# locus\tchromosome\tposition\n")
        for nm, (c, pos) in zip(names, coords):
            fh.write("%s\t%s\t%d\n" % (nm, c, pos))

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
        chrom, pos = coords[l]
        lines.append("\t".join([chrom, str(pos), names[l], "T", alt, ".",
                                 "PASS", ".", "GT"] + calls))
    text = "\n".join(lines) + "\n"

    with open(prefix + ".vcf", "w", newline="\n") as fh:
        fh.write(text)
    if gzip_vcf:
        with gzip.open(prefix + ".vcf.gz", "wt", newline="\n") as fh:
            fh.write(text)

    return {"dlines": ploidy * sum(nind), "loci": nloci,
            "nd_rows": 1, "nd_cols": 2, "sort_by": 2,
            "file": os.path.basename(prefix) + ".stru",
            "samples": os.path.basename(prefix) + ".samples",
            "map": os.path.basename(prefix) + ".map",
            "coords": coords, "names": names}


def write_phased(prefix, nloci, npops=4, nind=8, nall=2, missing=0.0, seed=1,
                 ploidy=2, haploid_samples=(), chroms=1, spacing=1000,
                 indels=False, short_gt=(), unphased=(), gzip_vcf=False):
    """Write a phased VCF whose haplotypes are known independently of ADZE.

    The haplotype matrix is generated first -- one sequence of allele indices
    per gene copy, along the whole chromosome -- and the VCF is written from
    it, so a checker can assemble the expected haplotype of any gene copy
    without parsing a single record the way ADZE does. The matrix is written
    beside the VCF as <prefix>.haps:

        SAMPLE  COPY  <locus> <locus> ...
        ind0_0  0     0       1       .

    with '.' for an uncalled copy. One uncalled call makes that copy's whole
    block uncalled, so these dots are what the block-level missing rule acts
    on, and a checker reads them rather than inferring them.

    Deliberate irregularities, each switched on per fixture:
      haploid_samples  (pop, individual) pairs carrying one gene copy
      short_gt         (record, pop, individual): GT written with one fewer
                       allele than the sample has, which makes the trailing
                       copy uncalled there -- the ploidy shortfall rule
      unphased         records written with '/' instead of '|', for the
                       refusal test
      indels           ALT alleles of differing length, so a haplotype label
                       cannot be read as a fixed-width string

    Returns the metadata the suites need, including the per-copy matrix.
    """
    rng = random.Random(seed)
    if isinstance(nind, int):
        nind = [nind] * npops
    assert len(nind) == npops

    coords = place_loci(nloci, chroms, spacing)
    names = ["S%d" % l for l in range(nloci)]

    # copies[(pop, individual)] = list of per-locus allele indices, one list
    # per gene copy. A population's own allele frequencies differ, so private
    # haplotypes exist rather than being a tie.
    copies = {}
    for p in range(npops):
        freq = [rng.random() * 0.6 + 0.2 for _ in range(nloci)]
        for i in range(nind[p]):
            n = 1 if (p, i) in haploid_samples else ploidy
            per_copy = []
            for _ in range(n):
                seq = []
                for l in range(nloci):
                    if rng.random() < missing:
                        seq.append(None)
                    elif nall <= 2:
                        seq.append(1 if rng.random() < freq[l] else 0)
                    else:
                        seq.append(rng.randrange(nall))
                per_copy.append(seq)
            copies[(p, i)] = per_copy

    samples = [(p, i) for p in range(npops) for i in range(nind[p])]

    # A shortened GT does not convey its trailing copy at that record, so the
    # matrix says uncalled there too: .haps is the record of what the VCF
    # conveys, not of what was drawn, and a checker that read the drawn value
    # would expect a haplotype the file never carried.
    for (l, p, i) in short_gt:
        seq = copies[(p, i)]
        if len(seq) > 1:
            seq[-1][l] = None

    with open(prefix + ".haps", "w", newline="\n") as fh:
        fh.write("SAMPLE\tCOPY\t" + "\t".join(names) + "\n")
        for (p, i) in samples:
            for c, seq in enumerate(copies[(p, i)]):
                fh.write("ind%d_%d\t%d\t%s\n" % (
                    p, i, c, "\t".join("." if a is None else str(a) for a in seq)))

    with open(prefix + ".samples", "w", newline="\n") as fh:
        fh.write("# sample\tgrouping\n")
        for (p, i) in samples:
            fh.write("ind%d_%d\tPOP%d\n" % (p, i, p))

    alts = ["AT" * (a + 1) if indels else "ACGT"[a % 4]
            for a in range(1, max(2, nall))]
    header = ["##fileformat=VCFv4.2",
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
              "#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                               "FILTER", "INFO", "FORMAT"]
                              + ["ind%d_%d" % (p, i) for (p, i) in samples])]

    lines = list(header)
    for l in range(nloci):
        sep = "/" if l in unphased else "|"
        calls = []
        for (p, i) in samples:
            seq = copies[(p, i)]
            drop = (l, p, i) in short_gt
            out = [("." if seq[c][l] is None else str(seq[c][l]))
                   for c in range(len(seq) - (1 if drop and len(seq) > 1 else 0))]
            calls.append(sep.join(out))
        chrom, pos = coords[l]
        lines.append("\t".join([chrom, str(pos), names[l], "G",
                                ",".join(alts[:max(1, nall - 1)]), ".",
                                "PASS", ".", "GT"] + calls))

    text = "\n".join(lines) + "\n"
    with open(prefix + ".vcf", "w", newline="\n") as fh:
        fh.write(text)
    if gzip_vcf:
        with gzip.open(prefix + ".vcf.gz", "wt", newline="\n") as fh:
            fh.write(text)

    return {"file": os.path.basename(prefix) + ".vcf",
            "samples": os.path.basename(prefix) + ".samples",
            "haps": os.path.basename(prefix) + ".haps",
            "names": names, "coords": coords, "nloci": nloci,
            "groups": ["POP%d" % p for p in range(npops)],
            "copies": {("ind%d_%d" % (p, i)): copies[(p, i)] for (p, i) in samples}}


def read_haps(path):
    """The per-copy haplotype matrix back from a .haps file.

    {(sample, copy): [allele index or None per locus]}, plus the locus names.
    """
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        names = head[2:]
        out = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[(f[0], int(f[1]))] = [None if a == "." else int(a) for a in f[2:]]
    return names, out


# Phased fixtures for the haplotype suites. Each one isolates something that
# block assembly has to get right.
PHASED = [
    ("hap_basic",   dict(npops=4, nind=6, nloci=24, nall=2, seed=11)),
    ("hap_missing", dict(npops=4, nind=6, nloci=24, nall=2, missing=0.10, seed=12)),
    ("hap_multi",   dict(npops=3, nind=5, nloci=18, nall=5, seed=13, indels=True)),
    ("hap_haploid", dict(npops=3, nind=6, nloci=18, nall=2, seed=14, ploidy=1)),
    ("hap_mixed",   dict(npops=3, nind=6, nloci=18, nall=2, seed=15,
                         haploid_samples=((0, 0), (1, 2)))),
    ("hap_short",   dict(npops=3, nind=5, nloci=12, nall=2, seed=16,
                         short_gt=((3, 0, 1), (7, 2, 4)))),
    ("hap_twochrom", dict(npops=3, nind=6, nloci=20, nall=2, seed=17, chroms=2)),
    ("hap_uneven",  dict(npops=4, nind=[3, 6, 9, 12], nloci=16, nall=3, seed=18)),
    ("hap_unphased", dict(npops=2, nind=4, nloci=8, nall=2, seed=19, unphased=(2,))),
]


def build_phased(outdir):
    os.makedirs(outdir, exist_ok=True)
    return {name: write_phased(os.path.join(outdir, name), **kw)
            for name, kw in PHASED}


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
    with open(path, "w", newline="\n") as fh:
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
