#!/usr/bin/env python3
"""Every translation unit declares what it uses.

A standard-library function is only guaranteed to be declared if the header
that defines it is included. Implementations differ in what they pull in
transitively: libc++ hands you <cstring> through several other headers, so a
file that calls memcmp() without including it compiles on macOS and fails on
a Linux box with libstdc++ --- which is exactly how ADZE_main_tools.cpp came
to use memcmp() with no #include <cstring>, building here and failing on a
cluster.

This suite reads the sources rather than compiling them, so it catches that
on any machine, including the one that would have compiled it anyway. For
each file it collects the headers visible to it --- its own includes, plus
those of the project headers it includes, transitively --- and requires one
for every standard symbol it uses.

    ./test/includes.py [--src src]

Exit status is 0 only if every file declares everything it uses.
"""
import argparse
import os
import re
import sys

# symbol -> the header that must be visible. Only symbols whose header is
# unambiguous are listed: a name that is also a plausible member function
# (find, remove, sort on a container) would produce false positives, so the
# patterns below require the free-function spelling.
NEEDS = [
    (r"(?<![\w.>])\b(memcmp|memcpy|memmove|memset|strlen|strcmp|strncmp|strcpy|strncpy|strchr|strstr|strtok)\s*\(",
     ("<cstring>", "<string.h>")),
    (r"(?<![\w.>])\b(atoi|atof|atol|atoll|strtol|strtoll|strtod|getenv|exit|malloc|calloc|realloc|free|qsort|abort)\s*\(",
     ("<cstdlib>", "<stdlib.h>")),
    (r"(?<![\w.>])\b(printf|fprintf|sprintf|snprintf|fopen|fclose|fread|fwrite|fflush|rename)\s*\(",
     ("<cstdio>", "<stdio.h>")),
    (r"(?<![\w.>])\bstd::remove\s*\(|(?<![\w.>])(?<!std::)\bremove\s*\(\s*\w+\.c_str\(\)",
     ("<cstdio>", "<stdio.h>")),
    (r"(?<![\w.>])\b(isspace|isdigit|isalpha|isalnum|toupper|tolower)\s*\(",
     ("<cctype>", "<ctype.h>")),
    (r"(?<![\w.>])\b(sqrt|pow|log|log10|exp|floor|ceil|fabs)\s*\(",
     ("<cmath>", "<math.h>")),
    (r"(?<![\w.>])\b(time|clock|difftime|localtime|strftime)\s*\(",
     ("<ctime>", "<time.h>")),
    (r"(?<![\w.>])\b(isatty|fileno|access|unlink)\s*\(",
     ("<unistd.h>", "<stdio.h>")),
    (r"\bnumeric_limits\s*<", ("<limits>",)),
    (r"\bstd::sort\s*\(|(?<![\w.>])\bsort\s*\(\s*\w+\.begin\(\)",
     ("<algorithm>",)),
    (r"\bostringstream\b|\bistringstream\b|\bstringstream\b", ("<sstream>",)),
    (r"\bifstream\b|\bofstream\b", ("<fstream>",)),
    (r"\bvector\s*<", ("<vector>",)),
    (r"\blist\s*<", ("<list>",)),
    (r"\bunordered_map\s*<", ("<unordered_map>",)),
    (r"(?<!unordered_)\bmap\s*<", ("<map>",)),
    (r"\bsetw\s*\(|\bsetprecision\s*\(", ("<iomanip>",)),
    (r"\bgzopen\s*\(|\bgzFile\b", ("<zlib.h>",)),
    (r"\bstat\s*\(\s*\w+|\bstruct stat\b", ("<sys/stat.h>",)),
]

COMMENT = re.compile(r"//[^\n]*|/\*.*?\*/", re.S)
STRING = re.compile(r'"(?:[^"\\]|\\.)*"')
INCLUDE = re.compile(r'^\s*#\s*include\s+([<"][^>"]+[>"])', re.M)


def includes(path):
    return set(INCLUDE.findall(open(path).read()))


def visible(name, src, seen=None):
    """Headers this file can rely on: its own, plus its project headers'."""
    seen = seen or set()
    if name in seen or not os.path.exists(os.path.join(src, name)):
        return set()
    seen.add(name)
    out = set()
    for inc in includes(os.path.join(src, name)):
        if inc.startswith('"'):
            out |= visible(inc.strip('"'), src, seen)
        else:
            out.add(inc)
    return out


def body(path):
    """File contents with comments and string literals removed."""
    text = COMMENT.sub(" ", open(path).read())
    return STRING.sub('""', text)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    names = sorted(f for f in os.listdir(args.src) if f.endswith((".cpp", ".h")))
    built = compiled(args.src, names)

    npass = 0
    failures = []
    for name in names:
        if name not in built:
            continue                      # not part of the program
        text = body(os.path.join(args.src, name))
        seen = visible(name, args.src)
        for pattern, headers in NEEDS:
            hit = re.search(pattern, text)
            if not hit:
                continue
            if any(h in seen for h in headers):
                npass += 1
                if args.verbose:
                    print("PASS  %-22s %-16s (%s)" % (name, headers[0], hit.group(0).strip()))
            else:
                failures.append((name, headers[0], hit.group(0).strip()))
                print("FAIL  %-22s uses %s but no %s is visible"
                      % (name, hit.group(0).strip(), " or ".join(headers)))

    print("\n%d passed, %d failed  (%d files)" % (npass, len(failures), len(built)))
    return 1 if failures else 0


def compiled(src, names):
    """The files the build actually compiles, from the Makefile's SOURCES.

    A retired file left in the tree is not held to this: it is not built, so
    a missing include in it cannot break anyone's build.
    """
    text = open(os.path.join(src, "Makefile")).read()
    m = re.search(r"^SOURCES\s*=\s*((?:[^\n\\]|\\\n)*)", text, re.M)
    if not m:
        return set(names)
    sources = set(m.group(1).replace("\\\n", " ").split())
    headers = {n for n in names if n.endswith(".h")}
    return sources | headers


if __name__ == "__main__":
    sys.exit(main())
