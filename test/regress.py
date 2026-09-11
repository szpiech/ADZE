#!/usr/bin/env python3
"""Differential regression test: candidate ADZE build vs. a reference build.

Every case is run twice, in separate scratch directories, and all output files
are compared byte-for-byte.  A case passes only if the two builds produce the
same file set with the same bytes (after masking the wall-clock lines that are
expected to differ between runs) and the same exit status.

One qualification since 2.0 sweeps from g = 1: ADZE 1.0 started at g = 2 and
cannot produce a g = 1 row, so whole files no longer match.  The comparison is
over the rows both builds compute -- every row at g >= 2, byte-for-byte -- and
the g = 1 rows are accounted for separately: each results file must carry
exactly one per label.  Report it that way rather than as a file-level
identity.

    ./test/regress.py --candidate src/adze --reference /path/to/adze-1.0

Cases where the reference build dies on a signal are not comparable: ADZE 1.0
aborts on some legal inputs (see review section 3.1).  Those are reported as
XFAIL while the candidate still shares the defect and as FIXED once it does not,
and they never fail the suite.

Exit status is 0 only if every comparable case passes.
"""
import argparse
import filecmp
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen_data  # noqa: E402

# Lines that legitimately differ between two runs of the same build.
TIME_LINE = re.compile(r"\(d:h:m:s\)")


def cases(datasets):
    """(name, dataset, paramfile kwargs) for every regression case."""
    out = []
    for ds in sorted(datasets):
        base = dict(g=6, tol=1, full_r=1, full_p=1, full_c=1)
        out.append(("%s.default" % ds, ds, base))
        out.append(("%s.maxg2" % ds, ds, dict(base, g=2)))
        out.append(("%s.maxg12" % ds, ds, dict(base, g=12)))
        out.append(("%s.tol_strict" % ds, ds, dict(base, tol=0)))
        out.append(("%s.tol_mid" % ds, ds, dict(base, tol=0.2)))
        out.append(("%s.summary_only" % ds, ds,
                    dict(base, full_r=0, full_p=0, full_c=0)))
        out.append(("%s.comb2" % ds, ds, dict(base, comb=1, k="2")))
        out.append(("%s.comb_range" % ds, ds, dict(base, comb=1, k="1-3")))
        out.append(("%s.comb_list" % ds, ds, dict(base, comb=1, k="1,3")))
    return out


def run_one(binary, workdir, dataset_file, meta, prefix, kwargs, extra_argv=()):
    os.makedirs(workdir, exist_ok=True)
    shutil.copy(dataset_file, workdir)
    par = os.path.join(workdir, "case.par")
    gen_data.write_paramfile(par, meta, os.path.basename(dataset_file), prefix, **kwargs)
    proc = subprocess.run([os.path.abspath(binary), "case.par", *extra_argv],
                          cwd=workdir, capture_output=True, text=True)
    return proc


def outputs(workdir, dataset_file):
    skip = {os.path.basename(dataset_file), "case.par"}
    return sorted(f for f in os.listdir(workdir) if f not in skip)


def stat_g(fields):
    """The g a results row reports, or None if the row has no g.

    The label ahead of g is one field for a grouping and k fields for a tuple,
    so g is the first field that parses as an integer. Grouping names in the
    generated fixtures and in the distributed example are never purely
    numeric, which is what makes that unambiguous here; it is a convenience
    for the harness, not a rule the program imposes on labels.
    """
    for i, f in enumerate(fields):
        try:
            return int(f)
        except ValueError:
            continue
    return None


def masked(path, drop_g1=False):
    """File contents with the parts that legitimately differ removed.

    Timing lines differ between any two runs. In the run summary, blank lines
    are dropped too, so that a build reporting extra phases (e.g. a total
    runtime line the reference never wrote) still compares on the parameters
    it echoes.

    drop_g1 removes the g = 1 rows, which ADZE 1.0 never computed: its sweeps
    started at g = 2. Passing it for the candidate only is what keeps this a
    comparison of the rows both builds produce, rather than a looser one. The
    g = 1 rows are accounted for separately, in g1_rows below.
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    try:
        text = raw.decode()
    except UnicodeDecodeError:
        return raw
    lines = [l for l in text.splitlines() if not TIME_LINE.search(l)]
    if path.endswith("_summary") or path.endswith(".summary.txt"):
        lines = [l for l in lines if l.strip()]
    elif drop_g1:
        lines = [l for l in lines if stat_g(l.split()) != 1]
    return "\n".join(lines).encode()


def g1_rows(path):
    """(g = 1 rows, distinct labels) in a results file.

    A statistics file carries one row per label per g, so a build sweeping
    from g = 1 must add exactly one row per label -- no more, and not none.
    """
    labels = set()
    n = 0
    for line in open(path):
        fields = line.split()
        g = stat_g(fields)
        if g is None:
            continue
        label = " ".join(fields[:fields.index(str(g))])
        labels.add(label)
        if g == 1:
            n += 1
    return n, labels


# ADZE 1.0 returned -1 (exit status 255) for every failure. The candidate
# distinguishes usage, I/O and data-validation errors, so a 255 from the
# reference matches any of those.
CANDIDATE_FAILURE_CODES = (2, 3, 4, 255)


def comparable_status(ref_rc, cnd_rc):
    if ref_rc == cnd_rc:
        return True
    return ref_rc == 255 and cnd_rc in CANDIDATE_FAILURE_CODES


def reference_has_no_loci(directory, files):
    """True if the reference reported statistics over zero loci.

    ADZE 1.0 carried on when filtering removed every locus and printed rows of
    'nan -0 nan' with NUM_LOCI = 0.  Those cases have no comparable output; the
    candidate is expected to diagnose the input instead.  Result rows are
    'group g num_loci ...', so the third field carries the locus count.
    """
    for f in files:
        if f.endswith("_summary") or f.endswith("_deletedloci"):
            continue
        try:
            with open(os.path.join(directory, f)) as fh:
                lines = fh.read().splitlines()
        except (OSError, UnicodeDecodeError):
            continue
        for line in lines:
            fields = line.split()
            if len(fields) >= 3 and fields[2] == "0" and "NUM_LOCI" not in line:
                return True
    return False


# Files that carry statistics rows, and so gain a g = 1 row in the candidate.
def is_results_file(name):
    return not (name.endswith("_summary") or name.endswith(".summary.txt")
                or name.endswith("deletedloci"))


def compare(ref_dir, cnd_dir, files):
    """Files that differ beyond what is expected, and missing g = 1 rows.

    The candidate sweeps from g = 1 and the reference from g = 2, so the
    comparison is over the g >= 2 rows; the extra rows are checked for
    separately, since a build that quietly stopped writing them would
    otherwise pass.
    """
    diffs, missing = [], []
    for f in files:
        a, b = os.path.join(ref_dir, f), os.path.join(cnd_dir, f)

        if masked(a) != masked(b, drop_g1=is_results_file(f)):
            diffs.append(f)
            continue

        if is_results_file(f):
            n, labels = g1_rows(b)
            if n != len(labels):
                missing.append("%s: %d g=1 rows for %d labels"
                               % (f, n, len(labels)))
    return diffs, missing


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate", required=True, help="binary under test")
    ap.add_argument("--reference", required=True, help="known-good binary (ADZE 1.0)")
    ap.add_argument("--workdir", default=os.path.join(HERE, "work"))
    ap.add_argument("--data", default=None, help="dataset directory (default: <workdir>/data)")
    ap.add_argument("--filter", default=None, help="only run cases matching this substring")
    ap.add_argument("--keep", action="store_true", help="keep scratch directories")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    data_dir = args.data or os.path.join(args.workdir, "data")
    if args.data is None:
        shutil.rmtree(args.workdir, ignore_errors=True)
    metas = gen_data.build_all(data_dir)

    # The shipped example is the one case with externally published expected
    # output, so it is always included.
    example = os.path.join(HERE, os.pardir, "example")
    all_cases = cases(metas)
    if os.path.isdir(example):
        metas["example"] = {"dlines": 22, "loci": 5, "nd_rows": 1, "nd_cols": 5,
                            "sort_by": 5, "file": "small_data.stru"}
        shutil.copy(os.path.join(example, "small_data.stru"), data_dir)
        for tag, kw in [("published", dict(g=8, tol=0.2, comb=1, k="2",
                                           full_r=1, full_p=1, full_c=1)),
                        ("nofilter", dict(g=8, tol=1, comb=1, k="1-3",
                                          full_r=1, full_p=1, full_c=1))]:
            all_cases.append(("example.%s" % tag, "example", kw))

    npass = nfail = nskip = 0
    failures = []
    for name, ds, kw in all_cases:
        if args.filter and args.filter not in name:
            continue
        meta = metas[ds]
        dfile = os.path.join(data_dir, meta["file"])
        ref_dir = os.path.join(args.workdir, name, "ref")
        cnd_dir = os.path.join(args.workdir, name, "cnd")
        shutil.rmtree(os.path.join(args.workdir, name), ignore_errors=True)
        ref = run_one(args.reference, ref_dir, dfile, meta, "out", kw)
        cnd = run_one(args.candidate, cnd_dir, dfile, meta, "out", kw)

        if ref.returncode < 0:
            # Reference died on a signal: nothing to compare against.
            if cnd.returncode < 0:
                status, detail = "XFAIL", (
                    "reference and candidate both die on signal %d "
                    "(pre-existing defect)" % -ref.returncode)
            elif cnd.returncode != 0:
                status, detail = "FAIL", (
                    "reference died on signal %d, candidate exited %d" % (
                        -ref.returncode, cnd.returncode))
            else:
                status, detail = "FIXED", (
                    "reference died on signal %d, candidate exited 0 with %d outputs" % (
                        -ref.returncode, len(outputs(cnd_dir, dfile))))
        elif reference_has_no_loci(ref_dir, outputs(ref_dir, dfile)):
            # Reference reported statistics over zero loci; only the
            # candidate's handling of the degenerate input is checkable.
            if cnd.returncode == 0:
                status, detail = "XFAIL", (
                    "reference reported 0 loci; candidate also exited 0")
            else:
                status, detail = "FIXED", (
                    "reference reported statistics over 0 loci, candidate "
                    "exited %d with a diagnostic" % cnd.returncode)
        elif not comparable_status(ref.returncode, cnd.returncode):
            status, detail = "FAIL", "exit status %d vs %d" % (ref.returncode, cnd.returncode)
        elif ref.returncode != 0:
            # Both refused the input. What a failed run leaves behind is not a
            # contract, so only the refusal itself is compared.
            status, detail = "PASS", "both rejected the input (%d / %d)" % (
                ref.returncode, cnd.returncode)
        else:
            ref_files, cnd_files = outputs(ref_dir, dfile), outputs(cnd_dir, dfile)
            if ref_files != cnd_files:
                status, detail = "FAIL", "output file set differs: %s vs %s" % (
                    ref_files, cnd_files)
            else:
                diffs, missing = compare(ref_dir, cnd_dir, ref_files)
                if diffs:
                    status, detail = "FAIL", "byte differences in %s" % ", ".join(diffs)
                elif missing:
                    status, detail = "FAIL", "g=1 rows wrong: %s" % "; ".join(missing)
                else:
                    status, detail = "PASS", (
                        "%d files identical over g>=2, g=1 rows present"
                        % len(ref_files))

        if status == "FAIL":
            nfail += 1
            failures.append((name, detail, ref_dir, cnd_dir))
        elif status in ("FIXED", "XFAIL"):
            nskip += 1
        else:
            npass += 1
        if args.verbose or status != "PASS":
            print("%-5s %-34s %s" % (status, name, detail))

    if not args.keep and not failures:
        shutil.rmtree(args.workdir, ignore_errors=True)

    print("\n%d passed, %d failed, %d not comparable (reference aborts)" % (
        npass, nfail, nskip))
    if failures:
        print("\nscratch kept for inspection:")
        for name, detail, ref_dir, cnd_dir in failures:
            print("  %s: %s\n    %s\n    %s" % (name, detail, ref_dir, cnd_dir))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
