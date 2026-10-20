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
  assembly   the haplotype each gene copy carries, the per-grouping counts
             they produce, and the order the allele slots come out in
  degenerate a block of one SNP is the SNP, so --hap-snps 1 must reproduce an
             ordinary run byte for byte
  labels     that two haplotypes sharing their digits stay two alleles
  reencode   the whole path, against a VCF in which each block has been
             rewritten as a single record carrying its haplotypes as alleles
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


def source_rows(path):
    """{block name: {grouping: (callable copies, [count per slot])}} + order."""
    out, order = {}, []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            name = f[0]
            if name not in out:
                out[name] = {}
                order.append(name)
            counts = [int(x) for x in f[7].split(",")] if f[7] else []
            out[name][f[3]] = (int(f[4]), counts)
    return out, order


def expected_block(haps, ranks, groups, block, ngroups, gindex):
    """What a block's allele table should hold, from the fixture's haplotypes.

    Returns the per-slot counts in ADZE's slot order: a slot is created the
    first time a haplotype is seen, scanning groupings in output order, then
    samples by their rank within the grouping, then gene copies.
    """
    seen, slots = {}, []
    for (sample, copy) in sorted(haps, key=lambda k: (gindex[groups[k[0]]],
                                                      ranks[k[0]], k[1])):
        seq = haps[(sample, copy)]
        if any(seq[i] is None for i in block):
            continue
        label = "|".join(str(seq[i]) for i in block)
        if label not in seen:
            seen[label] = len(slots)
            slots.append([0] * ngroups)
        slots[seen[label]][gindex[groups[sample]]] += 1
    return slots


def sample_ranks(vcf, groups, gorder):
    """Each sample's rank within its grouping, in VCF column order."""
    _, samples = records(vcf)
    seen = {g: 0 for g in gorder}
    out = {}
    for name in samples:
        g = groups.get(name)
        if g is None:
            continue
        out[name] = seen[g]
        seen[g] += 1
    return out


def case_assembly(s, adze, w):
    """Haplotypes, their counts, and the slot order they come out in."""
    for name, kw in gen_data.PHASED:
        if "unphased" in kw:
            continue
        meta = gen_data.write_phased(os.path.join(w, name), **kw)
        vcf = os.path.join(w, meta["file"])
        coords, _ = records(vcf)
        _, haps = gen_data.read_haps(os.path.join(w, meta["haps"]))
        groups = groups_of(os.path.join(w, meta["samples"]))
        gorder = meta["groups"]
        gindex = {g: i for i, g in enumerate(gorder)}
        ranks = sample_ranks(vcf, groups, gorder)

        for width in (3, 4):
            dump = os.path.join(w, "%s_a%d.src" % (name, width))
            r = run(adze, w, ["--data", meta["file"], "--samples", meta["samples"],
                              "--tolerance", "1", "--hap-snps", str(width),
                              "--dry-run"],
                    env={"ADZE_DUMP_SOURCE": dump})
            if r.returncode != 0:
                s.check("assembly.%s.%d.runs" % (name, width), False,
                        "exit %d: %s" % (r.returncode, r.stderr.strip()[:100]))
                continue

            rows, order = source_rows(dump)
            blocks = blocks_by_snps(coords, width)
            s.check("assembly.%s.%d.blocks" % (name, width),
                    len(order) == len(blocks),
                    "source yielded %d blocks, layout gives %d" % (len(order), len(blocks)))

            bad = []
            for b, block in enumerate(blocks):
                if b >= len(order):
                    break
                want = expected_block(haps, ranks, groups, block,
                                      len(gorder), gindex)
                got = rows[order[b]]
                for gi, g in enumerate(gorder):
                    mine = [row[gi] for row in want]
                    theirs = got.get(g, (0, []))[1]
                    if mine != theirs:
                        bad.append((order[b], g, mine, theirs))
            s.check("assembly.%s.snps%d" % (name, width), not bad,
                    "%d block-grouping tables differ, e.g. %s" % (len(bad), bad[:2]))


def outputs(workdir, prefix):
    """The result files a run wrote, excluding its parameter echo.

    The summary legitimately differs: it records the flags the run was given,
    and one of these runs was given --hap-snps.
    """
    return sorted(f for f in os.listdir(workdir)
                  if f.startswith(prefix + ".") and "summary" not in f)


def case_degenerate(s, adze, w):
    """--hap-snps 1 reproduces an ordinary run byte for byte.

    A block of one SNP is that SNP: the same gene copies are callable, the
    haplotype label is the allele itself, and the block takes the record's own
    name. So every output file must match an ordinary run exactly -- the
    strongest invariant available here, since it holds the whole pipeline to
    the one it replaces rather than to an expectation of mine.
    """
    runs = [
        ("plain",   []),
        ("pops",    ["--pops", "POP0,POP1"]),
        ("windows", ["--window-bp", "6000", "--step-bp", "3000",
                     "--min-window-loci", "2"]),
        ("tuples",  ["--combinations", "--tuples-k", "1,2", "--full-tuples"]),
        ("full",    ["--full-richness", "--full-private", "--max-g", "4"]),
        ("atg",     ["--at-g", "max"]),
        ("legacy",  ["--legacy", "--full-richness"]),
    ]

    for name, kw in gen_data.PHASED:
        if "unphased" in kw:
            continue
        meta = gen_data.write_phased(os.path.join(w, name), **kw)
        base = ["--data", meta["file"], "--samples", meta["samples"],
                "--tolerance", "1"]

        for tag, extra in runs:
            a = "%s_%s_snp" % (name, tag)
            b = "%s_%s_hap" % (name, tag)
            ra = run(adze, w, base + extra + ["--out-prefix", a])
            rb = run(adze, w, base + extra + ["--hap-snps", "1",
                                              "--out-prefix", b])
            s.check("degenerate.%s.%s.status" % (name, tag),
                    ra.returncode == rb.returncode,
                    "exit %d without, %d with --hap-snps 1" % (
                        ra.returncode, rb.returncode))
            if ra.returncode != 0:
                continue

            fa, fb = outputs(w, a), outputs(w, b)
            s.check("degenerate.%s.%s.files" % (name, tag),
                    [f[len(a):] for f in fa] == [f[len(b):] for f in fb],
                    "file sets differ: %s vs %s" % (fa, fb))

            diff = [f[len(a):] for f in fa
                    if open(os.path.join(w, f), "rb").read()
                    != open(os.path.join(w, b + f[len(a):]), "rb").read()]
            s.check("degenerate.%s.%s.bytes" % (name, tag), not diff,
                    "differs in %s" % diff)

    # The comparison would be worthless if it passed for any block width, so
    # one run with real blocks must differ from the SNP-level run.
    meta = gen_data.write_phased(os.path.join(w, "hap_basic"),
                                 **dict(gen_data.PHASED)["hap_basic"])
    base = ["--data", meta["file"], "--samples", meta["samples"], "--tolerance", "1"]
    run(adze, w, base + ["--out-prefix", "ctl_snp"])
    run(adze, w, base + ["--hap-snps", "2", "--out-prefix", "ctl_hap"])
    same = [f[len("ctl_snp"):] for f in outputs(w, "ctl_snp")
            if open(os.path.join(w, f), "rb").read()
            == open(os.path.join(w, "ctl_hap" + f[len("ctl_snp"):]), "rb").read()]
    s.check("degenerate.control", not same,
            "blocks of 2 SNPs produced identical %s, so the comparison proves "
            "nothing" % same)


def case_labels(s, adze, w):
    """Two haplotypes that share their digits must stay two alleles.

    A label is the copy's allele indices joined by '|'. Drop the separator and
    the haplotypes 1,11 and 11,1 both read as "111", merging two alleles into
    one -- which is visible in the counts, so this case is a real check rather
    than an inspection of label text (ADZE never prints a label).

    Written by hand: twelve alleles drawn at random almost never put a
    digit-sharing pair in the same block.
    """
    vcf = os.path.join(w, "labels.vcf")
    with open(vcf, "w", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                  "FILTER", "INFO", "FORMAT", "s1", "s2"]) + "\n")
        alt = ",".join("ACGT"[i % 4] * (i + 1) for i in range(11))
        # s1 copy 0 carries 1 then 11; s2 copy 0 carries 11 then 1.
        fh.write("\t".join(["chr1", "100", "r1", "G", alt, ".", "PASS", ".",
                            "GT", "1|1", "11|11"]) + "\n")
        fh.write("\t".join(["chr1", "200", "r2", "G", alt, ".", "PASS", ".",
                            "GT", "11|11", "1|1"]) + "\n")
    with open(os.path.join(w, "labels.samples"), "w", newline="\n") as fh:
        fh.write("s1\tP1\ns2\tP1\n")

    dump = os.path.join(w, "labels.src")
    r = run(adze, w, ["--data", "labels.vcf", "--samples", "labels.samples",
                      "--tolerance", "1", "--hap-snps", "2", "--dry-run"],
            env={"ADZE_DUMP_SOURCE": dump})
    s.check("labels.runs", r.returncode == 0, "exit %d" % r.returncode)
    if r.returncode != 0:
        return

    rows, order = source_rows(dump)
    counts = rows[order[0]]["P1"][1] if order else []
    s.check("labels.distinct", sorted(counts, reverse=True) == [2, 2],
            "1,11 and 11,1 should be two alleles of two copies each; got %s" % counts)


def record_fields(vcf):
    """[(chromosome, position, id)] in file order, and the sample names."""
    rows, samples = [], []
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]
                continue
            rows.append((f[0], int(f[1]), f[2]))
    return rows, samples


def reencode(src, dst, haps, samples, blocks, minsnps=1):
    """Rewrite a phased VCF with one record per haplotype block.

    Each block becomes a single record whose alleles are the distinct
    haplotypes its gene copies carry, so running ADZE on this file WITHOUT
    haplotype mode must reproduce the haplotype run on the original. The
    allele indices need not match ADZE's slot numbering: slots are created in
    first-seen order over the same samples in the same column order, so any
    consistent assignment gives the same ordering.

    A copy whose haplotype is uncalled is written '.', and every sample keeps
    its own number of gene copies, so the gene-copy totals and the
    missing-data denominators are the ones the original run had.
    """
    rows, _ = record_fields(src)
    ploidy = {}
    for (sample, copy) in haps:
        ploidy[sample] = max(ploidy.get(sample, 0), copy + 1)

    out = ["##fileformat=VCFv4.2",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
           "#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                            "FILTER", "INFO", "FORMAT"] + samples)]

    for block in blocks:
        if len(block) < minsnps:
            continue
        chrom, start, first_id = rows[block[0]]
        name = first_id if len(block) == 1 else "%s:%d-%d" % (
            chrom, start, rows[block[-1]][1])
        if len(block) == 1 and first_id == ".":
            name = "%s:%d" % (chrom, start)

        index, calls = {}, []
        for sample in samples:
            copies = []
            for c in range(ploidy.get(sample, 0)):
                seq = haps.get((sample, c))
                if seq is None or any(seq[i] is None for i in block):
                    copies.append(".")
                    continue
                label = tuple(seq[i] for i in block)
                if label not in index:
                    index[label] = len(index)
                copies.append(str(index[label]))
            calls.append("|".join(copies) if copies else ".")

        nalt = max(1, len(index) - 1)
        alt = ",".join("A" * (i + 1) for i in range(nalt))
        out.append("\t".join([chrom, str(start), name, "G", alt, ".", "PASS",
                              ".", "GT"] + calls))

    with open(dst, "w", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")
    return len(out) - 3


def case_reencode(s, adze, w):
    """The haplotype path against the ordinary path on re-encoded input."""
    shapes = [
        ("plain",   []),
        ("full",    ["--full-richness", "--full-private", "--max-g", "4"]),
        ("tuples",  ["--combinations", "--tuples-k", "1,2", "--full-tuples"]),
        ("pops",    ["--pops", "POP0,POP1"]),
        ("windows", ["--window-bp", "9000", "--step-bp", "9000",
                     "--min-window-loci", "1"]),
        ("tol",     ["--tolerance", "0.4"]),
    ]

    for name, kw in gen_data.PHASED:
        if "unphased" in kw:
            continue
        meta = gen_data.write_phased(os.path.join(w, name), **kw)
        vcf = os.path.join(w, meta["file"])
        coords, samples = records(vcf)
        _, haps = gen_data.read_haps(os.path.join(w, meta["haps"]))

        for width in (3, 5):
            blocks = blocks_by_snps(coords, width)
            flat = os.path.join(w, "%s_re%d.vcf" % (name, width))
            reencode(vcf, flat, haps, samples, blocks)

            for tag, extra in shapes:
                a = "%s_re%d_%s_hap" % (name, width, tag)
                b = "%s_re%d_%s_flat" % (name, width, tag)
                base = ["--samples", meta["samples"]]
                if "--tolerance" not in extra:
                    base += ["--tolerance", "1"]

                ra = run(adze, w, ["--data", meta["file"]] + base + extra +
                         ["--hap-snps", str(width), "--out-prefix", a])
                rb = run(adze, w, ["--data", os.path.basename(flat)] + base +
                         extra + ["--out-prefix", b])

                s.check("reencode.%s.%d.%s.status" % (name, width, tag),
                        ra.returncode == rb.returncode,
                        "haplotype run exited %d, re-encoded run %d" % (
                            ra.returncode, rb.returncode))
                if ra.returncode != 0:
                    continue

                fa, fb = outputs(w, a), outputs(w, b)
                s.check("reencode.%s.%d.%s.files" % (name, width, tag),
                        [f[len(a):] for f in fa] == [f[len(b):] for f in fb],
                        "file sets differ: %s vs %s" % (fa, fb))

                diff = []
                for f in fa:
                    other = os.path.join(w, b + f[len(a):])
                    if not os.path.exists(other):
                        continue
                    if open(os.path.join(w, f), "rb").read() != open(other, "rb").read():
                        diff.append(f[len(a):])
                s.check("reencode.%s.%d.%s.bytes" % (name, width, tag), not diff,
                        "differs in %s" % diff)


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
    for case in (case_layout, case_callable, case_assembly,
                 case_degenerate, case_labels, case_reencode,
                 case_filter):
        case(s, args.candidate, args.workdir)

    print("\n%d passed, %d failed" % (s.passed, len(s.failed)))
    return 1 if s.failed else 0


if __name__ == "__main__":
    sys.exit(main())
