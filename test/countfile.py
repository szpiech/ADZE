#!/usr/bin/env python3
"""The converted counts of a STRUCTURE file: where they go, and who owns them.

A STRUCTURE file is converted to per-locus counts before the sweep. Without
--count-file the conversion is the run's own temporary file, beside the
output, removed when the run ends. With --count-file it is kept and reused.

  concurrent  two runs started together in one directory with the same
              --out-prefix, on different data, must each give the result
              they give alone. Before each run had its own temporary file,
              the second rewrote the first's counts while the first read
              them, and the first to finish removed the other's.
  cleanup     no temporary file is left after a run, whether it succeeds or
              fails after the conversion.
  named       a --count-file is kept, reused by the next run on the same
              data, converted again for different data, and two runs that
              convert into the same named file at once both give the right
              result and leave no partial file behind.

The concurrent cases depend on timing, so each is run at several delays
between the two starts. Against the build before each run had its own
temporary file, the cases at the three shorter delays fail (one of them by
exiting 0 with the wrong numbers), and at the longest the first run finishes
before the second begins.

    ./test/countfile.py --candidate src/adze
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_data

BASE = ["--group-col", "2", "--tolerance", "1", "--max-g", "10", "--quiet"]
DELAYS = (0.0, 0.05, 0.2, 0.5)


class Suite:
    def __init__(self, verbose):
        self.verbose, self.passed, self.failed = verbose, 0, []

    def check(self, name, ok, detail=""):
        if ok:
            self.passed += 1
            if self.verbose:
                print("PASS  %s" % name)
        else:
            self.failed.append(name)
            print("FAIL  %s  %s" % (name, detail))


def start(adze, w, args):
    return subprocess.Popen([os.path.abspath(adze)] + args, cwd=w,
                            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                            text=True)


def run(adze, w, args):
    return subprocess.run([os.path.abspath(adze)] + args, cwd=w,
                          capture_output=True, text=True)


def same(a, b):
    with open(a, "rb") as fa, open(b, "rb") as fb:
        return fa.read() == fb.read()


def temporaries(w):
    return sorted(glob.glob(os.path.join(w, "*.tmp")))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--workdir", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "work-countfile"))
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    w = args.workdir
    shutil.rmtree(w, ignore_errors=True)
    os.makedirs(w)
    # Large enough that a conversion is still running when the second run
    # starts; different shapes, so a mixed-up result cannot look right.
    gen_data.write_dataset(os.path.join(w, "a.stru"), npops=8, nind=20,
                           nloci=8000, nall=15, missing=0.02, seed=71)
    gen_data.write_dataset(os.path.join(w, "b.stru"), npops=5, nind=14,
                           nloci=6000, nall=7, missing=0.03, seed=72)

    s = Suite(args.verbose)
    A = args.candidate

    # Each data set alone, as the reference.
    for d in ("a", "b"):
        r = run(A, w, ["--data", d + ".stru"] + BASE + ["--out-prefix", "solo_" + d])
        s.check("solo.%s" % d, r.returncode == 0, r.stderr.strip()[-200:])

    # concurrent: same --out-prefix, different data
    for i, delay in enumerate(DELAYS):
        pa = start(A, w, ["--data", "a.stru"] + BASE + ["--out-prefix", "same",
                          "--out-richness", "ra%d" % i, "--out-private", "pa%d" % i])
        time.sleep(delay)
        pb = start(A, w, ["--data", "b.stru"] + BASE + ["--out-prefix", "same",
                          "--out-richness", "rb%d" % i, "--out-private", "pb%d" % i])
        ea, eb = pa.communicate()[1], pb.communicate()[1]
        for tag, p, e, mine, solo in (("a", pa, ea, "ra%d" % i, "solo_a"),
                                      ("b", pb, eb, "rb%d" % i, "solo_b")):
            ok = p.returncode == 0 and same(os.path.join(w, mine),
                                            os.path.join(w, solo + ".richness"))
            s.check("concurrent.delay%g.%s" % (delay, tag), ok,
                    "exit %d %s" % (p.returncode, e.strip()[-160:]))
        s.check("concurrent.delay%g.cleanup" % delay, not temporaries(w),
                "left behind: %s" % temporaries(w))

    # cleanup after a failure that comes after the conversion
    with open(os.path.join(w, "bad.tuples"), "w") as fh:
        fh.write("NO_SUCH_GROUPING\n")
    r = run(A, w, ["--data", "a.stru"] + BASE + ["--tuples", "bad.tuples",
                   "--out-prefix", "fail"])
    s.check("cleanup.after_failure", r.returncode != 0 and not temporaries(w),
            "exit %d, left behind: %s" % (r.returncode, temporaries(w)))

    # named: kept, reused, converted again for other data
    r = run(A, w, ["--data", "a.stru"] + BASE + ["--count-file", "k.cnt",
                   "--out-prefix", "k1"])
    s.check("named.kept", r.returncode == 0 and os.path.exists(os.path.join(w, "k.cnt")))
    r = subprocess.run([os.path.abspath(A), "--data", "a.stru", "--group-col", "2",
                        "--tolerance", "1", "--max-g", "10", "--count-file", "k.cnt",
                        "--out-prefix", "k2"], cwd=w, capture_output=True, text=True)
    s.check("named.reused", "Reusing the counts in k.cnt" in r.stderr
            and same(os.path.join(w, "k2.richness"), os.path.join(w, "solo_a.richness")),
            r.stderr.strip()[-200:])
    r = run(A, w, ["--data", "b.stru"] + BASE + ["--count-file", "k.cnt",
                   "--out-prefix", "k3"])
    s.check("named.reconverted_for_other_data", r.returncode == 0 and
            same(os.path.join(w, "k3.richness"), os.path.join(w, "solo_b.richness")))

    # named: two runs converting into the same file at once
    os.remove(os.path.join(w, "k.cnt"))
    for i, delay in enumerate(DELAYS):
        if os.path.exists(os.path.join(w, "k.cnt")):
            os.remove(os.path.join(w, "k.cnt"))
        p1 = start(A, w, ["--data", "a.stru"] + BASE + ["--count-file", "k.cnt",
                          "--out-prefix", "n1_%d" % i])
        time.sleep(delay)
        p2 = start(A, w, ["--data", "a.stru"] + BASE + ["--count-file", "k.cnt",
                          "--out-prefix", "n2_%d" % i])
        e1, e2 = p1.communicate()[1], p2.communicate()[1]
        ok = (p1.returncode == 0 and p2.returncode == 0 and
              same(os.path.join(w, "n1_%d.richness" % i), os.path.join(w, "solo_a.richness")) and
              same(os.path.join(w, "n2_%d.richness" % i), os.path.join(w, "solo_a.richness")))
        s.check("named.concurrent.delay%g" % delay, ok,
                "exit %d/%d %s" % (p1.returncode, p2.returncode, (e1 + e2).strip()[-160:]))
        s.check("named.concurrent.delay%g.cleanup" % delay, not temporaries(w),
                "left behind: %s" % temporaries(w))

    print("\n%d passed, %d failed" % (s.passed, len(s.failed)))
    if not s.failed:
        shutil.rmtree(w, ignore_errors=True)
    return 1 if s.failed else 0


if __name__ == "__main__":
    sys.exit(main())
