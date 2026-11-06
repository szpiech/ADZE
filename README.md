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


## License

ADZE is free software, distributed under the GNU General Public License,
version 3 or any later version. See `LICENSE` for the full text. Copyright
(C) 2008-2026 Zachary A Szpiech.

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
