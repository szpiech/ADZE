# Test suite

Seven suites. Four are driven by `gen_data.py`. `regress.py` holds the numbers to
ADZE 1.0; `equivalence.py` holds a build to another 2.x build, byte for byte;
`formats.py` holds the input formats to each other; `windows.py`
holds the windowed statistics to the layout they claim to cover. `docs.py`
needs no data and no binary: it holds the README to the program's own option
table. All run under `ctest` when the tree is configured with CMake
(`regression` needs `-DADZE_REFERENCE=/path/to/adze-1.0`, `equivalence` needs
`-DADZE_EQUIV_REFERENCE=/path/to/adze-2.x`).

## ADZE_DUMP_COUNTS

Not a suite, but the facility the reader work is checked with. Setting
`ADZE_DUMP_COUNTS=<file>` makes a run write every count the reader produced —
one row per locus and grouping, with that grouping's total, its missing gene
copies and the per-allele counts in slot order — before any locus is filtered:

    ADZE_DUMP_COUNTS=counts.tsv ./src/adze --data d.vcf --samples d.samples ...

Two builds that produce the same dump read the same data, so a change to how
the file is read can be checked where the difference would be, instead of
waiting for it to surface as a difference in a mean. It is an environment
variable rather than a flag because it is a development facility and not part
of the interface. The STRUCTURE and VCF readers produce identical dumps for
the same genotypes, coordinates included, which is `formats.py`'s guarantee
stated one level lower.

## equivalence.py

For a change that is meant to alter nothing a user can see — a rewrite of how
the data is read or swept, rather than of what is computed. It runs the
candidate and a 2.x reference over the same inputs with the same flags and
requires every output file to be byte-identical, masking only the wall-clock
lines that differ between any two runs:

    ./test/equivalence.py --candidate src/adze --reference /path/to/adze-pre

119 cases: seven fixture shapes (including uneven sample sizes, haploid data,
40 alleles per locus, and a two-chromosome set for windows) crossed with the
invocations that exercise each output path — the ladder and `--at-g`, both
layouts, `--full-*`, tuples by k, a statistic subset, a grouping filter,
windows by locus and by basepair, all three input encodings, and `--dry-run`.
Exit status and stdout are compared too.

A deliberate change is declared rather than tolerated, and the declaration is
still an exact comparison: `--fulldata-locus-major` rewrites the *reference's*
`_fulldata` into the layout decision 1 of `doc/streaming.md` calls for and
requires the candidate to equal it byte for byte, so a changed value still
fails. The suite is negative-controlled — a perturbed statistic, a deleted
output file, shuffled `_fulldata` rows, the declared layout offered without
the flag, and an altered value under the flag are each caught.

## regress.py

A differential test: it runs a candidate build and a reference build over the
same inputs and requires the outputs to agree. It is how the 2.0 performance
and interface work is kept numerically faithful to ADZE 1.0.

"Agree" is exact wherever the two builds are supposed to produce the same
bytes, and explicitly not where they are not:

| part of the output | compared |
|---|---|
| labels, g, locus counts | exactly |
| per-locus values in `_fulldata` | exactly, cell by cell — the layouts are transposes |
| MEAN, VAR, STD_ERR | `1e-9` relative, values below `1e-12` treated as zero |
| `_deletedloci` | byte-for-byte |
| rows at g = 1 | not compared: 1.0 has none. Counted instead — one per label per file |

The summary tolerance exists because 2.0 accumulates the mean and variance in
one pass where 1.0 summed twice over the stored per-locus values: 1.0 prints a
variance of `2.95823e-31` on the distributed example where every deviation is
zero, and 2.0 prints `0`. `1e-9` is four orders tighter than the six
significant digits that are printed, so any visible disagreement still fails —
checked by perturbing the sixth digit of a mean and of a variance, both caught.

Build a reference binary once from the 1.0 sources (tag/commit `a70f765`), then:

    make -C src                                  # or: cmake --build build
    ./test/regress.py --candidate src/adze --reference /path/to/adze-1.0

83 cases run: nine synthetic datasets (`gen_data.py`) crossed with nine
parameter settings, plus the shipped `example/` dataset in both its published
configuration and an unfiltered one. Cases where the reference build aborts on a
signal are not comparable and are reported as XFAIL (candidate shares the
defect) or FIXED (candidate does not).

## formats.py

VCF input has no 1.0 counterpart to compare against, so it is checked against
the STRUCTURE reader instead. The same genotypes are written three ways --
`.stru`, `.vcf` and `.vcf.gz` -- and the three runs must be byte-identical:

    ./test/formats.py --candidate src/adze

That comparison is exact rather than approximate because both readers order a
locus's allele slots by the same rule, so every statistic sums its terms in the
same order. Five dataset shapes are covered, including unequal sample sizes and
haploid data. The compressed case is skipped, with a message, when the build
has no zlib.

## windows.py

Sliding windows have no reference implementation either, so each case is an
invariant:

    ./test/windows.py --candidate src/adze

A single window over a single-chromosome dataset must reproduce the
genome-wide output exactly; every window's mean, variance and standard error
must match a recomputation from the per-locus `*_fulldata` values (read down
that g's column); windows
must tile as asked, stay inside one chromosome, and hold each interior locus
twice when the step is half the width; evenly spaced loci must give the same
statistics whether the window is set in basepairs or in loci; VCF and
STRUCTURE-plus-map must agree; thread count must not matter; and unsorted
positions, a chromosome split across blocks, a locus with no coordinate, and
STRUCTURE without a map must all be refused.

The recomputation compares means to 2e-5 and variances to 1e-3, because it
reads the six-significant-digit output rather than the doubles behind it, and
a variance built from small deviations amplifies that rounding. Both bounds
are far tighter than any real error would be.

## gen_data.py

Also usable standalone. `write_dataset` writes one STRUCTURE file;
`write_pair` writes the same genotypes as STRUCTURE, VCF and compressed VCF,
plus the sample map and the locus map, which is what `formats.py` and
`windows.py` consume. `place_loci` deals loci out over chromosomes at a fixed
spacing, so a window of k spacings holds exactly k loci and the expected
layout can be written down by hand.

    ./test/gen_data.py /tmp/adze-data

## haplotypes.py

Holds the haplotype-block machinery to haplotypes assembled from the
fixtures' own `.haps` matrices, which `gen_data.write_phased` writes before it
writes a single record — so nothing here is checked against ADZE's own account
of the haplotypes. Needs no reference binary. The cases, in rising order of
strength:

- *layout*, *callable*, *assembly*: block boundaries, spans, the gene copies
  each block can score under the all-or-nothing missing rule, and the
  per-grouping allele tables including slot order.
- *degenerate*: `--hap-snps 1` must reproduce an ordinary run byte for byte,
  across eight fixtures and seven run shapes. Controlled — blocks of two SNPs
  must differ, so the case cannot pass vacuously.
- *reencode*: each block is rewritten as a single VCF record carrying its
  haplotypes as alleles, and an ordinary run on that file must match the
  haplotype run byte for byte. This is the one that exercises the whole path
  against an independent implementation.
- *blocks*, *dryrun*, *closed*, *refusals*, *filter*: the `_blocks` report,
  the dry run's layout summary, the exact closed forms at one gene copy, the
  error paths, and which blocks `--tolerance` drops.

    ./test/haplotypes.py --candidate src/adze

## pool_test.cpp

Tests the thread pool in `src/ADZE_threads.h` on its own, with 1, 2, 3, 8 and
64 threads. Every item must be run exactly once, the worker index must be in
range, a sum computed by the pattern the sweep uses (per-item values, then a
fold in item order) must give the same bits at every pool size, and an
exception thrown by one item must reach the caller. CMake builds it as
`pool_test` and runs it under ctest as `pool`.

## threads.py

Checks that the number of threads does not change the output. Each case is
run with 1, 2, 3, 4, 8 and 64 threads, and every file the run writes must be
byte-identical to the file written with one thread. The only exception is the
run summary, from which the elapsed times and the thread count are removed
before the comparison. The cases cover VCF, compressed VCF, STRUCTURE and
phased VCF input, and every kind of output, including the window files, the
`_fulldata` files, the combinations and the `_blocks` file. Some runs set
`ADZE_SWEEP_BATCH` or `ADZE_TUPLE_BUDGET` to split the work into small
batches, so that batch boundaries fall inside windows and between passes over
the combinations. No reference binary is needed.

Two controls check that the comparison can fail. Every case must write at
least one non-empty result file, and runs with `MAX_G` of 5 and 6 must be
reported as different.

    ./test/threads.py --candidate src/adze

## includes.py

Reads the sources and requires each compiled file to include a header for
every standard-library symbol it uses --- following its project headers one
level, since a header it includes may supply one. It exists because
implementations differ in what they hand over transitively: `libc++` supplies
`<cstring>` through other headers, so a `memcmp()` with no
`#include <cstring>` built on macOS and failed on a Linux cluster with
`libstdc++`. Reading the source catches that anywhere, including on the
machine where it would have compiled.

    ./test/includes.py

## docs.py

The only one that needs neither a binary nor data. It parses the option table
out of `src/ADZE_pfile.cpp` and the option table out of `README.md` and fails
if they disagree: a README row naming an option the program no longer has, or
claiming a default the program does not apply. The comparison is permissive
about wording and strict about the value, and it is one-directional on purpose
-- accepting the README's text inside the table's default instead would pass a
README saying `1` where the program uses `0.1`. Confirmed by running it against
a deliberately wrong README.

## In CI

`.github/workflows/ci.yml` runs all four on pushes to `master` and `devel` and
on pull requests into either, and
builds the 1.0 reference from `master` in the same run rather than carrying a
binary around, so `regress.py` always has something to compare against. The
suites job runs them through `ctest` on Linux and macOS; the Makefile job runs
them directly, so it needs nothing but `make` and `python3`. On Windows
`regress.py` does not run: 1.0 built with mingw-w64 reports no statistics at
all — its `putNj` stores its value and then falls off the end of a `bool`
function, and in that build the store does not survive — so there is no usable
reference there. The Windows job runs the other three and re-runs the
distributed example against its committed output instead. OpenMP is present
on the Linux runner and absent on the macOS one, which is how the
thread-invariance case gets exercised in one place and reports itself skipped
in the other.
