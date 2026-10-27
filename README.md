# ADZE

[![CI](https://github.com/szpiech/ADZE/actions/workflows/ci.yml/badge.svg)](https://github.com/szpiech/ADZE/actions/workflows/ci.yml)

Allelic Diversity Analyzer

ADZE is a program that implements the rarefaction method for analyzing allelic diversity across populations while correcting for sample size differences. Using individual multilocus genotype data on genetic polymorphisms, ADZE computes estimates of allelic richness, private allelic richness, and private allelic richness for combinations of populations. 

ZA Szpiech, M Jakobsson, NA Rosenberg. (2008) ADZE: a rarefaction approach for counting alleles private to combinations of populations. Bioinformatics 24: 2498-2504. 

## Binaries

Tagged releases include pre-built archives, which are built and checked by
`.github/workflows/release.yml` from the tagged commit:

| archive | contents |
|---|---|
| `adze-<version>-linux-x86_64.tar.gz` | statically linked, with zlib |
| `adze-<version>-macos-universal2.tar.gz` | x86\_64 and arm64 in one binary, with zlib |
| `adze-<version>-windows-x86_64.zip` | statically linked, with zlib |

Each archive contains the binary, this README, the manual and the `example/`
directory, and nothing else needs to be installed. The Linux and Windows
binaries are statically linked, and the macOS binary only loads libraries in
`/usr/lib`. `SHA256SUMS` covers all three archives.

Before an archive is published, the binary in it is used to re-run the
distributed example, the output is compared with the committed output, and the
format suite is run against the same binary.

## Building

Only a C++17 compiler is required. If zlib is available, it is used to read
compressed input, and the program builds and runs without it. `--threads` uses
threads from the C++ standard library, so it needs no other library.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cmake --install build --prefix ~/.local     # optional
```

Or with `make`:

```sh
make -C src                        # no compressed input
make -C src ZLIB=1                 # with compressed input
```

`EXTRA_CXXFLAGS` and `EXTRA_LDFLAGS` add flags for a particular build, such as
a static link or several architectures, and they are used to build the release
archives:

```sh
make -C src ZLIB=1 EXTRA_LDFLAGS=-static
make -C src ZLIB=1 EXTRA_CXXFLAGS="-arch x86_64 -arch arm64" \
                   EXTRA_LDFLAGS="-arch x86_64 -arch arm64"
```

Note that setting `CXXFLAGS` on the command line replaces the value in the
Makefile, and `make` then ignores the `+=` lines that `ZLIB=1` depends on. For example, `ZLIB=1 CXXFLAGS=...` builds without support for
compressed input, and no warning is given. The `EXTRA_` variables avoid this
problem.

On Windows, the same commands work in the MINGW64 shell of MSYS2. Install
`mingw-w64-x86_64-gcc`, and also `mingw-w64-x86_64-cmake` and
`mingw-w64-x86_64-ninja` to build with CMake. The compiler adds `.exe` to an
output name that does not have it, so the `make` build needs
`make -C src TARGET=adze.exe`. Input files with Windows line endings are read
correctly, and the format suite checks this on every platform.

## Running

```sh
adze --data mydata.stru --group-col 5 --out-prefix mypops
adze --data cohort.vcf.gz --samples populations.tsv --out-prefix cohort
adze --data cohort.vcf.gz --samples populations.tsv \
     --window-bp 100000 --step-bp 25000 --out-prefix scan
adze --help
```

### Input formats

ADZE reads two input formats, chosen with `--format` or, by default, from the
file name:

- The STRUCTURE format used by version 1.0, with one row per gene copy, a
  number of label columns, and then one allele code per locus.
- VCF, uncompressed or compressed with `gzip` or `bgzip`. Each record is one
  locus, and each allele index is one allele type. A sample contributes as many
  gene copies as its `GT` field holds, so haploid and diploid records can be in
  the same file, and only `GT` is read. Because a VCF does not give population
  labels, `--samples FILE` is required. This file has two whitespace-separated
  columns, the sample name and then the grouping name. A `#` starts a comment,
  and any further columns are ignored.

The same genotypes in either format give byte-identical results, and
`test/formats.py` checks this.

### Sliding windows

Any of the three statistics can be reported in windows along the genome, in
addition to once for the whole data set. The size of a window is given either
in basepairs (`--window-bp`) or in loci (`--window-loci`). Basepair windows
start at position 1, so the windows do not depend on where the loci are.
Locus windows count the loci that remain after filtering. By default, each
`--step-*` is equal to the window size, so the windows do not overlap. No
window spans two chromosomes.

Each statistic then writes an additional `*_windows` file, which always has a
header and is always tab-separated, with the columns `CHROM START END
POP_GROUPING G NUM_LOCI MEAN VAR STD_ERR`. The files for the whole data set do
not change.

Every locus must have a position. A VCF gives the positions, and STRUCTURE
input needs `--loci-map`. A locus without a position is an error. `MAX_G` is
determined from the whole genome, so every window covers the same range of g
and the windows can be compared with each other. The loci must be in genome
order, and if they are not, ADZE reports the first locus that is out of order
rather than sorting the file.

### Standardized sample sizes

Every statistic is reported at each g from 1 to `MAX_G`. At g = 1, allelic
richness is exactly 1, because one gene copy always carries one allele, so this
row is a check on how the alleles were counted rather than a measurement.
Private allelic richness at g = 1 is the probability that a single gene copy
drawn from the grouping carries an allele that is missed by a single draw from
each of the other groupings.

The statistics files are tab-separated with a header row, and a statistic that
is undefined at a given g is written as `NA`. `--legacy` writes the version
1.0 layout instead, with no header, space-separated columns, and undefined
rows left out rather than marked. Leaving out these rows is the reason the
default changed. A grouping that is undefined at every g does not appear in the
version 1.0 layout at all, so a run can finish and report nothing for that
grouping without any explanation. The window files are always tab-separated
with a header, and the `_deletedloci` and `_summary` files are not affected by
`--legacy`.

The `_fulldata` files are transposed relative to version 1.0, with one row for
each locus and grouping and one column for each g. They are tab-separated with
a header, and an undefined value is written as `NA`. This is the same in both
layouts, so `--legacy` does not reproduce these files. In the version 1.0
layout (one row for each g, one column for each locus), the file cannot be
written until every locus has been computed, and at genome scale it has more
columns than most tools can read. The locus is in the first column because the
values for every grouping at a locus are computed together, before the next
locus is read. The "Locus-specific output" section of the manual shows both
layouts.

`_deletedloci` lists the dropped loci in the order of the input (version 1.0
listed them in descending order of position), and it is written even when
every locus is dropped.

The default `--tolerance` is 0.1, while version 1.0 used 1. If every locus is
kept, a single locus at which one grouping has no scored gene copies reduces
`MAX_G` to 1 and makes that grouping undefined at every g, so the run finishes
without reporting any statistics. `--tolerance 1` gives the version 1.0
behavior, and a parameter file that sets `TOLERANCE` is not affected. When a
grouping cannot be used, ADZE reports it by name, with the number of loci
involved, and `--dry-run` reports the same information before computing any
statistics.

`--at-g N` reports a single g, and `--at-g max` reports g = `MAX_G` without
the value needing to be known in advance. If the request is larger than
`MAX_G`, or larger than the smallest number of gene copies scored at any locus,
ADZE uses that value instead and prints a warning. Most of the cost of a full
run is in accumulating and writing the statistics at every g, so reporting one
g is much faster. For example, on 20 000 loci in 12 groupings with all 66 pairs
of groupings and `MAX_G` of 33, the full run took 3.42 s, `--at-g 2` took
0.60 s and `--at-g 33` took 0.65 s.

`LOCI`, `DATA_LINES` and `NON_DATA_COLS` are measured from the input file, and
by default `MAX_G` is the largest standardized sample size that the loci
remaining after filtering support. A run therefore only needs a data file and,
for STRUCTURE input with several label columns, the label column that names the
grouping. `--dry-run` reports the measured dimensions, the sample size of each
grouping, the largest feasible `MAX_G` and the number of combinations an
analysis would evaluate, and then exits without computing any statistics.

Parameter files from version 1.0 still work, with the same keywords and the
same short flags, and options given on the command line take precedence over
the file:

```sh
adze small_paramfile.txt -t 0.2
```

Messages, warnings and progress bars are written to standard error, and the
results are written to the output files. The exit status is 0 on success, 2
for a usage error, 3 for an input or output error, and 4 when the data file is
inconsistent or does not match the parameters.

### Options

The defaults listed here are the values the program uses. `adze --help` prints
the same values, since both are read from the table of options in
`src/ADZE_pfile.cpp`.

| option | effect | default |
|---|---|---|
| `--data FILE` | genotype file: STRUCTURE layout or VCF, optionally gzipped | required |
| `--out-prefix P` | write `P.richness`, `P.private`, `P.tuples_k*`, `P.summary.txt` | `adze` |
| `--max-g N` | largest standardized sample size | the largest the data supports |
| `--tolerance X` | drop loci where any grouping exceeds fraction X missing | `0.1` (`1` keeps all, as in 1.0) |
| `--group-col N` | which label column names the grouping | the last label column |
| `--missing STR` | code for a missing allele | `-9` |
| `--loci N`, `--data-lines N`, `--non-data-cols N` | declare a dimension instead of measuring it | detected |
| `--non-data-rows N` | header rows before the genotypes | `1` |
| `--stat richness,private,tuples` | compute only the named statistics | `richness,private` |
| `--pops A,B,C` | analyse only these groupings | every grouping |
| `--exclude-pops X` | analyse everything except these groupings | none excluded |
| `--combinations` | also compute private alleles of grouping tuples | off |
| `--tuples-k 1-3` | tuple sizes to enumerate (`-k`, with `--combinations`) | none |
| `--tuples FILE` | private alleles of the tuples named in FILE, one per line | every k-subset |
| `--at-g N` / `--at-g max` | report a single g instead of every g from 1 to `MAX_G` | every g |
| `--samples FILE` | sample-to-grouping map, required for VCF input | none |
| `--format auto\|structure\|vcf` | input format; `auto` reads the file name | `auto` |
| `--loci-map FILE` | locus coordinates, required for windowed STRUCTURE input | none |
| `--hap-snps N` | collapse runs of N consecutive SNPs into haplotype alleles (phased VCF only) | off |
| `--hap-bp N` | collapse each N-basepair run of SNPs into haplotype alleles (phased VCF only) | off |
| `--min-hap-snps N` | drop a haplotype block holding fewer than N SNPs | 1 |
| `--count-file FILE` | keep the counts a STRUCTURE file is converted into, and reuse them next run | temporary, removed |
| `--convert-chunk N` | loci per conversion pass | as many as fit 64 MB |
| `--window-bp N` / `--window-loci N` | report each statistic in sliding windows, sized in basepairs or in loci | off |
| `--step-bp N` / `--step-loci N` | how far a window advances | the window width |
| `--min-window-loci N` | skip windows holding fewer than N loci | `1` |
| `--full-richness`, `--full-private`, `--full-tuples` | also write the per-locus values | off |
| `--legacy` | write version 1.0's output layout | off |
| `--threads N` | number of threads; the output does not depend on it | `1` |
| `--progress` | progress bars | on when stderr is a terminal |
| `--quiet` | suppress progress and informational messages | off |
| `--dry-run` | report the detected layout, groupings and feasible `MAX_G`, then stop | off |
| `--write-template [FILE]` | generate a commented parameter file | — |

### Continuous integration

`.github/workflows/ci.yml` runs on pushes to `master` and `devel`, and on any
pull request into either branch, so a topic branch is checked when it is
proposed rather than on every push:

| job | what it checks |
|---|---|
| `suites` | the test suites under `ctest`, on Linux and macOS, against a version 1.0 binary built from `master` in the same run |
| `plain make` | the two Makefile configurations (bare and `ZLIB=1`) build and pass the suites |
| `MSYS2 MINGW64` | the program builds and runs on Windows through MSYS2, and reproduces the committed output of the distributed example |
| `distributed example` | `adze small_paramfile.txt` in `example/` reproduces the tracked result files |
| `manual` | the manual builds from source, and its generated option reference agrees with the option table |
| `races` | the thread pool and the program, built with ThreadSanitizer, show no data race on the determinism suite |

The version 1.0 reference is built from `master` in every run rather than kept
as a binary, so the comparison with version 1.0 is repeated on every run. This
is the only step in CI that needs GSL.

The comparison with version 1.0 is only run on Linux and macOS. In version
1.0, `Population::putNj` stores its value and then reaches the end of a `bool`
function without returning, and the version 1.0 build hides the compiler
warning with `-w`. With mingw-w64, and with GCC 15 or later, the stored value
is lost, so every sample size is 0, every statistic is undefined, and version
1.0 writes headers without numbers while reporting success. Before building
the reference, CI therefore runs `test/patch_ref10.py`, which adds the missing
`return` to `putNj` and changes nothing else. The patched reference gives
byte-identical output to an unpatched build in which the store is kept. The
Windows job does not yet use it, and instead re-runs the distributed example
and compares it with output committed from another platform.

## Memory

The genotypes are not held in memory. A run reads the input twice. The first
pass counts the gene copies for each grouping at each locus, which is enough to
determine the missing-data filter, the largest g for each grouping, `MAX_G`
and the windows. The second pass computes every statistic from one table of
*Q* per locus, and it keeps a running mean and variance rather than the
per-locus values. `--dry-run` only uses the first pass.

The following peak memory use was measured on simulated diploid data with 10
groupings of 20 individuals, both statistics and `--max-g 10` (median of three
runs):

| loci | VCF before | VCF now | STRUCTURE before | STRUCTURE now |
|---|---|---|---|---|
| 20 000 | 13 MB | 7 MB | 15 MB | 16 MB |
| 80 000 | 49 MB | 12 MB | 55 MB | 35 MB |
| 320 000 | 195 MB | 32 MB | 216 MB | 133 MB |
| 1 280 000 | 764 MB | 145 MB | 886 MB | 260 MB |

"Before" is the version 2.0 build from before this change, which held the data
in memory. VCF runs take the same time as before. STRUCTURE runs take 1.3 to
1.9 times as long, because the file has one row per gene copy and is first
converted to counts arranged by locus. The converted file is written in the
same directory as the output and removed when the run ends. It is small (an
80 000-locus file of 122 MB becomes 7.2 MB), and the block of loci held in
memory during the conversion is why the peak for STRUCTURE input is higher than
for VCF. `--count-file FILE` keeps the converted file, and the next run over
the same data reuses it. The file is not used, and the reason is reported, if
the data file has changed in size or date, or if the grouping filters, the
missing-data code or the column layout differ. `--convert-chunk N` sets the
number of loci converted in each pass.

## Haplotype alleles

With a phased VCF, `--hap-snps N` (or `--hap-bp N`) collapses consecutive SNPs
into blocks, and the haplotypes carried by the gene copies are counted instead
of the alleles at single SNPs:

```sh
adze --data cohort.vcf.gz --samples pops.tsv --hap-snps 5 --out-prefix hap
```

The estimators do not change, since a block is treated as a locus whose
alleles are haplotypes. However, this has three consequences:

- Phase and sort order are required. Every `GT` with more than one allele must
  use `|`, and a `/` stops the run with an error that gives the record and the
  sample. Records must be sorted by chromosome and position (`bcftools sort`),
  since a block is a run of consecutive records. The STRUCTURE format does not
  say which chromosome copy an allele is on, so it is not supported in this
  mode.
- One uncalled genotype makes the whole block uncalled for that gene copy, so
  missing data accumulate along a block. With 2% missing genotypes, a block of
  50 SNPs keeps about a third of its gene copies, which reduces the largest
  feasible `MAX_G`.
- Allelic richness saturates. In a long block, every gene copy carries a
  different haplotype, and allelic richness is then equal to g regardless of
  the data. Private allelic richness reaches the same limit. The run reports
  how many blocks are in this state, and the `_blocks` file (block, region,
  SNPs, distinct haplotypes, and callable gene copies in each grouping) shows
  how close each block is to it.

Haplotype blocks can be used with sliding windows. Each window averages over
the blocks whose first SNP it contains, and `--window-loci` counts blocks. A
block can extend past the end of its basepair window, and the run reports how
many do. Blocks always fall inside windows when `--hap-bp` divides both
`--window-bp` and `--step-bp`.

`--hap-snps 1` gives exactly the same output as an ordinary run, and the test
suite checks this. `--dry-run` reports the blocks that a given block size
produces before any statistic is computed.

## Documentation

`ADZE_Manual.pdf` describes the estimators, the input formats, every option,
the output files, and what changed from version 1.0. Its source files are in
`doc/`. The option reference is generated from the table of options in the
source code, so it always matches the options the program accepts:

```sh
make -C doc            # build doc/ADZE_Manual.pdf (needs pdflatex)
make -C doc install    # copy it to ADZE_Manual.pdf at the top of the repository
```

## Testing

`test/regress.py` compares the output of a build with the output of a
reference build on the same inputs. The labels, the numbers of loci and the
per-locus values must agree exactly, and the mean, variance and standard error
must agree to a relative tolerance of `1e-9`. Version 2.0 computes the mean and
variance in a single pass, so the last digits of a variance can differ (see
`test/README.md`). To run it, build a binary from the version 1.0 sources, and
then:

```sh
./test/regress.py --candidate src/adze --reference /path/to/adze-1.0
```

`test/haplotypes.py` checks the haplotype blocks against haplotypes assembled
from the simulated data without using ADZE. `test/bench.py` measures run time
and memory against a version 1.0 build (locus filtering, the sweep over g,
the combinations of groupings, and reading the data). See `test/README.md` for
the other suites.
