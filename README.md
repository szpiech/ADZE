# ADZE
Allelic Diversity Analyzer

ADZE is a program that implements the rarefaction method for analyzing allelic diversity across populations while correcting for sample size differences. Using individual multilocus genotype data on genetic polymorphisms, ADZE computes estimates of allelic richness, private allelic richness, and private allelic richness for combinations of populations. 

ZA Szpiech, M Jakobsson, NA Rosenberg. (2008) ADZE: a rarefaction approach for counting alleles private to combinations of populations. Bioinformatics 24: 2498-2504. 

## Building

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

| option | effect |
|---|---|
| `--out-prefix P` | write `P.richness`, `P.private`, `P.tuples_k*`, `P.summary.txt` |
| `--stat richness,private,tuples` | compute only the named statistics |
| `--pops A,B,C` / `--exclude-pops X` | restrict the analysis to some groupings |
| `--tuples FILE` | private alleles of the tuples named in FILE, one per line, instead of every k-subset |
| `--tuples-k 1-3` | tuple sizes to enumerate (`-k`, with `--combinations`) |
| `--samples FILE` | sample-to-grouping map, required for VCF input |
| `--format auto\|structure\|vcf` | input format; `auto` reads the file name |
| `--window-bp N` / `--window-loci N` | report each statistic in sliding windows, sized in basepairs or in loci |
| `--step-bp N` / `--step-loci N` | how far a window advances (default: its own width, i.e. no overlap) |
| `--min-window-loci N` | skip windows holding fewer than N loci |
| `--loci-map FILE` | locus coordinates, required for windowed STRUCTURE input |
| `--tolerance X` | drop loci where any grouping exceeds fraction X missing |
| `--tsv` | tab-separated output with a header row and `NA` for undefined values |
| `--threads N` | parallelize the per-locus loops (OpenMP builds) |
| `--quiet`, `--progress` | control what reaches stderr |
| `--write-template [FILE]` | generate a commented parameter file |

Results do not depend on `--threads`: loci are processed independently and
reduced in locus order, so any thread count gives bit-identical output.

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
both binaries over the same inputs and requires every output file to match
byte-for-byte. Build a binary from the 1.0 sources, then:

```sh
./test/regress.py --candidate src/adze --reference /path/to/adze-1.0
```

`test/bench.py` measures the same axes as the performance work (locus
filtering, the `MAX_G` sweep, tuple enumeration, read time and memory, thread
scaling). See `test/README.md`.
