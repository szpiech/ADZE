# ADZE

[![CI](https://github.com/szpiech/ADZE/actions/workflows/ci.yml/badge.svg)](https://github.com/szpiech/ADZE/actions/workflows/ci.yml)

Allelic Diversity Analyzer

ADZE is a program that implements the rarefaction method for analyzing allelic diversity across populations while correcting for sample size differences. Using individual multilocus genotype data on genetic polymorphisms, ADZE computes estimates of allelic richness, private allelic richness, and private allelic richness for combinations of populations. 

ZA Szpiech, M Jakobsson, NA Rosenberg. (2008) ADZE: a rarefaction approach for counting alleles private to combinations of populations. Bioinformatics 24: 2498-2504. 

## Binaries

Tagged releases carry pre-built archives, built and checked by
`.github/workflows/release.yml` on the tag itself:

| archive | contents |
|---|---|
| `adze-<version>-linux-x86_64.tar.gz` | statically linked, with zlib and OpenMP |
| `adze-<version>-macos-universal2.tar.gz` | x86\_64 and arm64 in one binary, with zlib |
| `adze-<version>-windows-x86_64.zip` | statically linked, with zlib and OpenMP |

Each holds the binary, this README, the manual and the `example/` directory,
and needs nothing installed: the Linux and Windows binaries are static, and
the macOS one loads only libraries in `/usr/lib`. `SHA256SUMS` covers all
three.

The macOS build is serial. Linking OpenMP there would tie the binary to a
Homebrew `libomp` that a machine downloading it will not have, so `--threads`
on macOS wants a local build — the third `make` line below.

Every archive is checked before it is published, with the binary that will
ship: it re-runs the distributed example and compares against the committed
output, and runs the format suite against itself.

## Building

On Windows, MSYS2's MINGW64 shell works with the same commands; install
`mingw-w64-x86_64-gcc`, and `mingw-w64-x86_64-cmake` with
`mingw-w64-x86_64-ninja` for the CMake route. The compiler appends `.exe` to an
output name that has none, so the Makefile route wants
`make -C src TARGET=adze.exe`. Input files with Windows line endings are read
as they are; the carriage return never reaches an allele label, which the
format suite checks on every platform.

Nothing is required beyond a C++17 compiler. OpenMP (for `--threads`) and zlib
(for compressed input) are used when the toolchain provides them, and the
program builds and runs without either.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cmake --install build --prefix ~/.local     # optional
```

Or with plain make:

```sh
make -C src                        # serial, no compressed input
make -C src OPENMP=1 ZLIB=1        # GCC or Linux clang
make -C src OPENMP=1 ZLIB=1 OMPFLAGS="-Xpreprocessor -fopenmp" OMPLIBS=-lomp
                                   # Apple clang, with libomp installed
```

`EXTRA_CXXFLAGS` and `EXTRA_LDFLAGS` add flags that belong to a build rather
than to the program — a static link, or several architectures — and are how
the release archives are built:

```sh
make -C src ZLIB=1 OPENMP=1 EXTRA_LDFLAGS=-static
make -C src ZLIB=1 EXTRA_CXXFLAGS="-arch x86_64 -arch arm64" \
                   EXTRA_LDFLAGS="-arch x86_64 -arch arm64"
```

Setting `CXXFLAGS` itself on the command line would replace the value in the
Makefile, and make then ignores the `+=` lines that `ZLIB=1` and `OPENMP=1`
rely on — so `ZLIB=1 CXXFLAGS=...` builds without compressed-input support and
says nothing about it. The `EXTRA_` hooks exist to avoid that.

## Running

```sh
adze --data mydata.stru --group-col 5 --out-prefix mypops
adze --data cohort.vcf.gz --samples populations.tsv --out-prefix cohort
adze --data cohort.vcf.gz --samples populations.tsv \
     --window-bp 100000 --step-bp 25000 --out-prefix scan
adze --help
```

### Input formats

Two are read, chosen by `--format` or, by default, by the file name:

- the STRUCTURE-like layout ADZE has always taken — one row per gene copy,
  label columns then one allele code per locus;
- VCF, plain or `gzip`/`bgzip` compressed. Each record is one locus and each
  allele index one allele type; a sample contributes as many gene copies as its
  `GT` field holds, so haploid and diploid records can mix, and only `GT` is
  read. Because a VCF carries no population labels, `--samples FILE` is
  required: two whitespace-separated columns, sample name then grouping name,
  `#` for comments, further columns ignored.

The same genotypes in either format give byte-identical results; `test/formats.py`
checks that.

### Sliding windows

Any of the three statistics can be reported along the genome instead of once
for the whole dataset. A window is measured either in basepairs (`--window-bp`,
anchored at position 1 so the intervals do not depend on where the loci fall)
or in surviving loci (`--window-loci`), never both; each `--step-*` defaults to
its window width, so the default is a non-overlapping tiling. No window crosses
a chromosome.

Each statistic then writes an extra `*_windows` file — always with a header,
always tab-separated — with columns `CHROM START END POP_GROUPING G NUM_LOCI
MEAN VAR STD_ERR`. The whole-dataset files are unchanged.

Windows need coordinates for every locus: VCF supplies them, STRUCTURE input
needs `--loci-map`, and a locus with no coordinate is an error rather than a
locus quietly left out. `MAX_G` stays resolved genome-wide so every window
covers the same range of *g* and the scan is comparable end to end. Loci must
be in genome order; ADZE says which locus breaks the order rather than sorting
the file for you.

### Standardized sample sizes

Every statistic is reported at each g from 1 to `MAX_G`. At g = 1 allelic
richness is exactly 1 by construction (one gene copy carries one allele), which
makes that row a check on the allele binning rather than a measurement; private
allelic richness at g = 1 is the chance that a single copy drawn from the
grouping carries an allele that single draws from every other grouping all
miss.

Result files are tab-separated with a header row, and a statistic undefined at
a given g is written as `NA`. `--legacy` is the only way to change that: it
asks for 1.0's layout — no
header, space separators, undefined rows omitted rather than marked. That
omission is why the default changed: a grouping undefined at every g disappears
from the legacy file entirely, so a run can complete and report nothing with no
visible reason. Window files were always tab-separated with a header, and `_deletedloci` and
`_summary` are unaffected either way. The `_fulldata` files are **transposed with respect to 1.0** — one row per
grouping and locus with a column per g, tab-separated with a header, `NA` where
the value is undefined — in both layouts, so `--legacy` does not reproduce that
file. 1.0's orientation (a row per g, a column per locus) cannot be written
until the whole sweep is finished and is unreadable at genome scale; the
manual's "Locus-specific output" section shows both and how to transpose back.

The `--tolerance` default is 0.1, not 1.0's `1`. Keeping every locus means a
single locus where one grouping scored nothing drives `MAX_G` to 1 and leaves
that grouping undefined at every g, so the run completes and reports nothing;
`--tolerance 1` restores the old behaviour, and a parameter file that declares
`TOLERANCE` is unaffected. When a grouping is left unusable, adze now names it
with the locus count rather than leaving you to infer it from empty files, and
`--dry-run` reports the same ceiling before computing anything.

`--at-g N` reports one g instead, and `--at-g max` reports at whatever `MAX_G`
resolves to without needing to know the number. A request above the reachable
ceiling — `MAX_G`, or the smallest number of gene copies scored anywhere — is
met at that ceiling with a warning. Since the recurrence in g is sequential the
sweep still climbs to the requested g, so the saving scales with how low it is:
on 20k loci in 12 groupings with all 66 pairwise tuples and `MAX_G` 33,
`--at-g 2` took 0.58s against 6.72s for the full ladder.

`LOCI`, `DATA_LINES` and `NON_DATA_COLS` are measured from the input file, and
`MAX_G` defaults to the largest standardized sample size the surviving loci
support, so the minimum invocation is a data file and (for multi-column
headers) which label column names the grouping. `--dry-run` reports the
detected layout, the per-grouping sample sizes, the largest feasible `MAX_G`
and the number of tuples an analysis would evaluate, without computing
anything.

Parameter files from ADZE 1.0 still work, with the same keywords and the same
short flags, and command-line options override the file:

```sh
adze small_paramfile.txt -t 0.2
```

Diagnostics, warnings and progress bars go to stderr; results go to the output
files. Exit status is 0 on success, 2 for a usage error, 3 for an I/O error and
4 when the data file disagrees with itself or with the parameters.

### Options worth knowing

Defaults are those the program actually applies; `adze --help` prints the same
values, both being read from the one option table in `src/ADZE_pfile.cpp`.

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
| `--window-bp N` / `--window-loci N` | report each statistic in sliding windows, sized in basepairs or in loci | off |
| `--step-bp N` / `--step-loci N` | how far a window advances | the window width |
| `--min-window-loci N` | skip windows holding fewer than N loci | `1` |
| `--full-richness`, `--full-private`, `--full-tuples` | also write the per-locus values | off |
| `--legacy` | write version 1.0's output layout | off |
| `--threads N` | parallelize the per-locus loops (OpenMP builds) | `1` |
| `--progress` | progress bars | on when stderr is a terminal |
| `--quiet` | suppress progress and informational messages | off |
| `--dry-run` | report the detected layout, groupings and feasible `MAX_G`, then stop | off |
| `--write-template [FILE]` | generate a commented parameter file | — |

### Continuous integration

`.github/workflows/ci.yml` runs on pushes to `master` and `devel`, and on any
pull request proposing to merge into either — so a topic branch is checked when
it is proposed rather than on every push:

| job | what it establishes |
|---|---|
| `suites` | all four suites under `ctest`, on Linux and macOS, against a version 1.0 binary built from `master` in the same run |
| `plain make` | the three Makefile configurations (bare, `ZLIB=1`, `ZLIB=1 OPENMP=1`) build and pass the suites |
| `MSYS2 MINGW64` | builds and runs on Windows through MSYS2, and reproduces the distributed example's committed output |
| `distributed example` | `adze small_paramfile.txt` in `example/` reproduces the tracked result files |
| `manual` | the manual builds from source and its generated option reference agrees with the option table |

The 1.0 reference is rebuilt from `master` on every run rather than kept as a
binary, so the byte-identity claim is re-established rather than inherited. It
is the only thing in CI that needs GSL.

The differential comparison runs on Linux and macOS only. ADZE 1.0's
`Population::putNj` stores its value and then falls off the end of a `bool`
function, which its own build hides with `-w`; in a mingw-w64 build that store
does not survive, so every sample size stays 0, every statistic is undefined,
and 1.0 writes headers without numbers while reporting success. There is
nothing to compare against there, so Windows instead re-runs the distributed
example and compares against output committed from another platform, which
pins this build's numbers without involving 1.0.

Results do not depend on `--threads`: loci are processed independently and
accumulated in locus order, so any thread count gives bit-identical output —
checked at 1, 2, 4 and 8 threads over the statistics, `_fulldata` and window
files.

## Documentation

`ADZE_Manual.pdf` documents the estimators, the input format, every option, the
output files, and what changed from version 1.0. Its sources are in `doc/`; the
option reference is generated from the option table in the source, so it cannot
drift from the parser:

```sh
make -C doc            # build doc/ADZE_Manual.pdf (needs pdflatex)
make -C doc install    # refresh ADZE_Manual.pdf at the repository root
```

## Testing

`test/regress.py` is a differential test against a reference build: it runs
both binaries over the same inputs and requires the outputs to agree — labels,
locus counts and per-locus values exactly, the three summary columns to `1e-9`
relative (2.0 accumulates them in one pass, so the last bits of a variance can
differ; see `test/README.md`). Build a binary from the 1.0 sources, then:

```sh
./test/regress.py --candidate src/adze --reference /path/to/adze-1.0
```

`test/bench.py` measures the same axes as the performance work (locus
filtering, the `MAX_G` sweep, tuple enumeration, read time and memory, thread
scaling). See `test/README.md`.
