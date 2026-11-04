#!/usr/bin/env python3
# ADZE: Allelic Diversity Analyzer
# Copyright (C) 2026 Zachary A Szpiech
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Make the ADZE 1.0 sources safe to build with a current compiler.

In version 1.0, Population::putNj stores its value and then reaches the end of
a bool function without returning. That is undefined behaviour, and its caller
never reads the result, so older compilers built it as intended. GCC 15 at -O2
(Ubuntu 26.04), gcc 16 even at -O0 and mingw-w64 do not: the stored value is
lost, and the reference then reports no statistics while exiting normally.

This adds the missing `return 1;` after the store, which is the only change.
It is the one function g++ -O2 -Wreturn-type reports in the 1.0 sources. A
build of the patched sources gives byte-identical output to an unpatched
build that does not exploit the undefined behaviour (gcc 11 at -O0), so the
reference still computes what version 1.0 computes.

    ./test/patch_ref10.py ref10/src

Exits nonzero, changing nothing, unless the store is found exactly once.
"""
import os
import sys

STORE = "      Nj[locus] = num;"


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    path = os.path.join(sys.argv[1], "ADZE_pop.cpp")
    with open(path, "rb") as fh:
        text = fh.read().decode("latin-1")

    eol = "\r\n" if "\r\n" in text else "\n"
    lines = text.split(eol)
    hits = [i for i, l in enumerate(lines) if l == STORE]
    if len(hits) != 1:
        sys.exit("%s: expected the putNj store once, found it %d times; "
                 "not patched" % (path, len(hits)))

    i = hits[0]
    if lines[i + 1].strip() == "return 1;":
        print("%s: already patched" % path)
        return 0
    if lines[i + 1].strip() != "}":
        sys.exit("%s: putNj no longer ends after the store; not patched" % path)

    lines.insert(i + 1, "      return 1;")
    with open(path, "wb") as fh:
        fh.write(eol.join(lines).encode("latin-1"))
    print("%s:%d: added 'return 1;' to Population::putNj" % (path, i + 2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
