#!/usr/bin/env python3
"""Hold the README's option table to the defaults the program actually applies.

The interface is described once, in the OPTIONS table in src/ADZE_pfile.cpp,
which the parsers, --help, --write-template and the manual's option reference
all read. The README is written by hand, so it is the one place that can drift:
this checks it cannot drift silently.

Three things are checked, none of which needs a built binary:

  exists     every long option the README names is in the OPTIONS table, so a
             renamed or removed flag cannot linger in the README
  default    where the README declares a default, it agrees with the table's
  covered    where the table gives a default, no README row claims a
             different one

A README row may cover several flags (one row for the three declared
dimensions, say) and may add a parenthetical note after the value, so the
comparison strips backticks and a trailing parenthetical and then asks whether
the table's default appears in what the README says: permissive about wording,
strict about the value. The containment is one-directional on purpose. Allowing
the README's text to be found inside the table's default instead would accept a
README saying 1 where the program uses 0.1, since "1" is a substring of "0.1"
-- which is exactly the drift this is meant to catch, and did catch when the
check was tried against a deliberately wrong README.

    ./test/docs.py [--readme README.md] [--source src/ADZE_pfile.cpp]
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

ENTRY = re.compile(
    r"\{\s*(\w+)\s*,\s*"
    r'(0|"[^"]*")\s*,\s*'
    r'(0|"[^"]*")\s*,\s*'
    r'(0|"[^"]*")\s*,\s*'
    r"(OPT_\w+)\s*,\s*"
    r'(0|"[^"]*")\s*,\s*'
    r'(0|"[^"]*")\s*,\s*'
    r'"([^"]*)"\s*,\s*'
    r'"((?:[^"\\]|\\.)*)"\s*\}',
    re.S)

# Ways the README may express "this option has no default value".
NO_DEFAULT = ("required", "none", "-", "--", "\u2014")

# Answered in main() before any parameter is parsed, so they are deliberately
# not in the OPTIONS table: --help and --version must work with no data file
# and no valid parameters, and --write-template writes the table out rather
# than reading from it.
NOT_IN_TABLE = ("--help", "--version", "--write-template")


def unquote(tok):
    return None if tok == "0" else tok[1:-1]


def options(path):
    """{long option: default string or None} from the OPTIONS table."""
    with open(path) as fh:
        src = fh.read()
    start = src.index("const OptSpec OPTIONS[] = {")
    end = src.index("\n};", start)
    out = {}
    for m in ENTRY.finditer(src[start:end]):
        out[unquote(m.group(3))] = unquote(m.group(7))
    if not out:
        raise SystemExit("no OPTIONS entries parsed from %s" % path)
    return out


def readme_rows(path):
    """(flags, default cell) for each row of a markdown table with a default column."""
    rows = []
    header = None
    for line in open(path):
        if not line.startswith("|"):
            header = None
            continue
        # A cell may contain an escaped pipe (an option listing alternatives),
        # which must not be read as a column separator.
        cells = [c.strip().replace("\\|", "|")
                 for c in re.split(r"(?<!\\)\|", line.strip().strip("|"))]
        if header is None:
            header = [c.lower() for c in cells]
            continue
        if set("".join(cells)) <= set("-: "):
            continue                       # the |---|---| separator
        if "default" not in header:
            continue
        flags = re.findall(r"--[a-z][a-z-]*", cells[0])
        rows.append((flags, cells[header.index("default")]))
    return rows


def normalize(text):
    text = text.replace("`", "")
    text = re.sub(r"\s*\([^)]*\)\s*$", "", text)   # drop a trailing note
    return text.strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--readme", default=os.path.join(ROOT, "README.md"))
    ap.add_argument("--source", default=os.path.join(ROOT, "src", "ADZE_pfile.cpp"))
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    table = options(args.source)
    rows = readme_rows(args.readme)
    if not rows:
        print("FAIL  no option table with a default column found in %s" % args.readme)
        return 1

    npass = 0
    failures = []

    def check(name, ok, detail):
        nonlocal npass
        if ok:
            npass += 1
            if args.verbose:
                print("PASS  %s" % name)
        else:
            failures.append((name, detail))
            print("FAIL  %s: %s" % (name, detail))

    for flags, cell in rows:
        check("row.names_an_option", bool(flags),
              "no long option in %r" % cell)
        said = normalize(cell)
        for flag in flags:
            if flag in NOT_IN_TABLE:
                continue
            check("exists.%s" % flag, flag in table,
                  "not in the OPTIONS table")
            if flag not in table:
                continue

            want = table[flag]
            if want is None:
                check("default.%s" % flag, said.lower() in NO_DEFAULT,
                      "table gives no default, README says %r" % said)
            else:
                check("default.%s" % flag, want.lower() in said.lower(),
                      "table says %r, README says %r" % (want, said))

    print("\n%d passed, %d failed  (%d README rows, %d options in the table)"
          % (npass, len(failures), len(rows), len(table)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
