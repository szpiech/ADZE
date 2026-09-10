# Test suite

Two suites, both driven by `gen_data.py`. `regress.py` holds the numbers to
ADZE 1.0; `formats.py` holds the input formats to each other. Both run under
`ctest` when the tree is configured with CMake (`regress` needs
`-DADZE_REFERENCE=/path/to/adze-1.0`).

## regress.py

A differential test: it runs a candidate build and a reference
build over the same inputs and requires the outputs to be byte-identical. It is
how the 2.0 performance and interface work is kept numerically faithful to
ADZE 1.0.

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

## gen_data.py

Also usable standalone. `write_dataset` writes one STRUCTURE file;
`write_pair` writes the same genotypes as STRUCTURE, VCF and compressed VCF
plus the sample map, which is what `formats.py` consumes.

    ./test/gen_data.py /tmp/adze-data
