# ADZE
Allelic Diversity Analyzer

ADZE is a program that implements the rarefaction method for analyzing allelic diversity across populations while correcting for sample size differences. Using individual multilocus genotype data on genetic polymorphisms, ADZE computes estimates of allelic richness, private allelic richness, and private allelic richness for combinations of populations. 

ZA Szpiech, M Jakobsson, NA Rosenberg. (2008) ADZE: a rarefaction approach for counting alleles private to combinations of populations. Bioinformatics 24: 2498-2504. 

## Building

No external dependencies; a C++17 compiler is enough. OpenMP is used when the
toolchain provides it.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cmake --install build --prefix ~/.local     # optional
```

Or with plain make:

```sh
make -C src                 # serial
make -C src OPENMP=1        # GCC or Linux clang
make -C src OPENMP=1 OMPFLAGS="-Xpreprocessor -fopenmp" OMPLIBS=-lomp
                            # Apple clang, with libomp installed
```

## Running

```sh
adze --data mydata.stru --group-col 5 --out-prefix mypops
adze --help
```

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
| `--tolerance X` | drop loci where any grouping exceeds fraction X missing |
| `--tsv` | tab-separated output with a header row and `NA` for undefined values |
| `--threads N` | parallelize the per-locus loops (OpenMP builds) |
| `--quiet`, `--progress` | control what reaches stderr |
| `--write-template [FILE]` | generate a commented parameter file |

Results do not depend on `--threads`: loci are processed independently and
reduced in locus order, so any thread count gives bit-identical output.

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
