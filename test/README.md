# Test suite

Three suites, all driven by `gen_data.py`. `regress.py` holds the numbers to
ADZE 1.0; `formats.py` holds the input formats to each other; `windows.py`
holds the windowed statistics to the layout they claim to cover. All three run
under `ctest` when the tree is configured with CMake (`regress` needs
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

## windows.py

Sliding windows have no reference implementation either, so each case is an
invariant:

    ./test/windows.py --candidate src/adze

A single window over a single-chromosome dataset must reproduce the
genome-wide output exactly; every window's mean, variance and standard error
must match a recomputation from the per-locus `*_fulldata` columns; windows
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
