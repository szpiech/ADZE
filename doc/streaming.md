# The streaming rewrite: what it changes, and what it must not

ADZE reads the whole dataset into memory and then sweeps it. Since the
statistics are accumulated in one pass (no per-locus values are stored any
more), nothing in the arithmetic needs more than one locus at a time, so the
dataset need not be resident at all. This branch makes the program read a
locus, compute from it, and forget it.

This document is the contract the work is held to. It was written before any
of it, so that a difference found later is either on this list or a defect.

## What must not change

Everything a user can see, apart from the three decisions below: the
statistics and their values, the *g* ladder, which rows appear in the
statistics files and in what order, the window output, the deleted-loci
report, the option names and their defaults, the exit codes, and the terms on
which version 1.0 is reproduced.

This is enforced rather than intended. `test/equivalence.py` runs the
candidate and a build from this branch's starting point over the same inputs
with the same flags and requires every output file to be byte-identical — 119
cases, covering each output path in all three input encodings. A declared
decision below is expressed *in* that suite as an exact expectation, not as a
relaxation: the expected file is computed from the reference's own output, so
a changed value still fails.

At a lower level, `ADZE_DUMP_COUNTS` (see `test/README.md`) writes the counts
the reader produced, per locus and grouping. A new reader is checked against
the old one there, where a miscount is a wrong number in a named locus rather
than a drifting mean.

## Decision 1: `_fulldata` becomes locus-major, with `LOCUS` first

A streaming sweep produces a locus's values for every grouping at once, so
the natural row order is by locus. Keeping today's grouping-major order would
mean buffering a stream per grouping and concatenating at the end — that is,
holding what the rewrite exists to stop holding.

Today:

    POP_GROUPING  LOCUS  G1  G2  G3
    POP0          L0     ...
    POP0          L1     ...
    POP1          L0     ...

After:

    LOCUS  POP_GROUPING  G1  G2  G3
    L0     POP0          ...
    L0     POP1          ...
    L1     POP0          ...

The columns swap because the first column is the sort key, and the file joins
to a locus table on it. Everything else about the file is unchanged: the same
rows, the same values, the same `NA` for a value undefined at that *g*, the
same tab separation and header, and `TUPLE` in place of `POP_GROUPING` in the
tuple file. `--legacy` does not restore the old order; as with the transpose
before it, that flag governs the statistics files.

Enforced by `test/equivalence.py --fulldata-locus-major`.

## Decision 2: the input is read twice, and the summary says so

Which loci survive `--tolerance`, what the default `MAX_G` resolves to, where
the windows fall and which loci were dropped are all properties of the whole
dataset, and all are needed before the first statistic can be computed. So a
scan comes first, keeping no per-locus values, and the sweep follows.

- **VCF**: the file is read twice.
- **STRUCTURE**: the file is read once, during the conversion of decision 3,
  which computes everything the scan needs as it goes; the sweep then reads
  the converted file.

There is no conditional single-pass mode. It would apply only when `--max-g`
or `--at-g` is given *and* nothing is filtered *and* no windows are requested,
which is a narrow case bought at the price of a second path through the
program for every future change to defend.

The run summary reports both phases, as it reports the phases today.

## Decision 3: STRUCTURE is converted, not streamed

A STRUCTURE file is individual-major: one row per gene copy, one column per
locus. A block of loci is therefore a column slice, and streaming one would
mean re-reading the file once per block. Rather than teach the engine two
input shapes, the engine gets one — a locus-major source — and STRUCTURE
input is converted to it: read once, written as temporary per-locus counts,
which the sweep then consumes exactly as it consumes a VCF.

Consequences worth stating: allele-label interning (the linear scan over
labels, the first-seen ordering) lives in the converter and leaves the hot
path; the conversion needs temporary disk, counts rather than genotypes; and
a dataset analysed repeatedly can be converted once.

VCF input pays none of this.

### The converted file

**What it holds.** A header naming the groupings and their gene copies, then
one record per locus: the locus name, its allele count, and the counts
slot-major with stride *J*. Coordinates are not repeated, the scan having
them already. Counts rather than genotypes, so it is small: an 80 000-locus
STRUCTURE file of 122 MB converts to 7.2 MB.

**Where it goes, and when it disappears.** Beside the output, as
`<out-prefix>.counts.tmp`, not in a system temporary directory — it can be
large, and the output directory is the one the user chose for large files.
It is removed when the run ends, whether or not the run succeeded. A file
the user names explicitly is never removed.

**Converting once for many runs.** A named file is reused when it describes
this run's data, which means all of: the same source path, the same size,
the same modification time, the same `--pops`/`--exclude-pops`, the same
missing-data code, the same grouping column, and the same header-row count —
plus the same groupings, in the same order, with the same gene copies. Every
one of those changes what the counts would be.

A mismatch is reported with its reason and the data converted again. It is
never used with a warning: a stale count file is yesterday's genotypes in
today's table, and no warning makes that acceptable.

**What the check cannot see.** Size and modification time are not a hash. An
edit that preserves both — rewriting a file in place within the filesystem's
timestamp resolution — would go unnoticed. Name a converted file when the
input is settled; while data is still moving, let the run convert into its
temporary file each time.

### Decision 4: the deleted-loci list is in file order

The sweep names a dropped locus as it goes past, because the scan does not
keep locus names -- a million of them would put back a per-locus cost this
rewrite exists to remove. The list is therefore in file order, where before
it was in descending index order. Same header sentence, same names, reversed.

A run where every locus is dropped stops before the sweep, and that is the
run whose user most needs the list, so it is written there too: from the
header row for STRUCTURE input, and by re-reading the records for a VCF,
which is affordable on a path that is about to exit.

### Decision 5: one pass, one completion line

Each statistic is still announced -- "Calculating allelic richness...", and
so on for private richness and each tuple size -- because each is still
being calculated. They are now calculated together, from one Q table per
locus, so there is one completion line rather than three, and the run
summary carries one "Statistics completed at" line in place of one per
statistic. No per-statistic timing is lost that meant anything: in a fused
pass there is no moment at which richness is finished and private richness
is not.
