# Test suite

`regress.py` is a differential test: it runs a candidate build and a reference
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

`gen_data.py` also serves as a standalone generator:

    ./test/gen_data.py /tmp/adze-data
