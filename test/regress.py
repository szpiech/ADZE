#!/usr/bin/env python3
"""Differential regression test: candidate ADZE build vs. a reference build.

Every case is run twice, in separate scratch directories, and the outputs are
compared: byte-for-byte for the deleted-loci list, cell by cell for the
transposed _fulldata files (the two layouts are transposes -- see the manual),
and field by field for the statistics files, where the label and locus-count
columns must match exactly and the three summary columns to a tolerance well
inside the printed precision (see SUMMARY_RTOL).  A case passes only if the two builds produce the
same file set with the same bytes (after masking the wall-clock lines that are
expected to differ between runs) and the same exit status.

Two qualifications. 2.0's default output layout is tab-separated with a header
and NA for undefined rows, so the candidate is run with --legacy to ask for the
layout 1.0 wrote; the flag is the only argv difference between the two runs.

And since 2.0 sweeps from g = 1: ADZE 1.0 started at g = 2 and
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

    # Kept beside the outputs: when a run misbehaves on a machine you cannot
    # reach, what it said about itself is usually the answer, and the scratch
    # directory is what CI uploads.
    with open(os.path.join(workdir, "console.out"), "w", newline="\n") as fh:
        fh.write(proc.stdout)
    with open(os.path.join(workdir, "console.err"), "w", newline="\n") as fh:
        fh.write(proc.stderr)

    return proc


def outputs(workdir, dataset_file):
    """Result files the run wrote, excluding the inputs and captured console."""
    skip = {os.path.basename(dataset_file), "case.par",
            "console.out", "console.err"}
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


# ADZE 1.0 returned -1 for every failure -- 255 as an exit code, and
# 4294967295 unsigned on Windows. The candidate distinguishes usage, I/O and
# data-validation errors, so that one reference status matches any of those.
CANDIDATE_FAILURE_CODES = (2, 3, 4, 255)


def comparable_status(ref_rc, cnd_rc):
    if ref_rc == cnd_rc:
        return True
    return (normal_status(ref_rc) == 255
            and normal_status(cnd_rc) in CANDIDATE_FAILURE_CODES)


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
def normal_status(returncode):
    """A process status as an exit code, with Windows' unsigned image undone.

    A program that calls exit(-1) -- which is how 1.0 reports every ordinary
    failure -- exits 255 on POSIX and 4294967295 on Windows, the same value
    seen as unsigned. Anything in the top 256 of the unsigned range is that,
    not a crash.
    """
    if returncode >= 0xFFFFFF00:
        return (returncode - 0x100000000) & 0xFF
    return returncode


def died(returncode):
    """True when the process was killed rather than exiting on its own.

    POSIX reports a signal as a negative status. Windows has no signals: a
    crash surfaces as the NTSTATUS value itself, and the failure severity
    occupies 0xC0000000 upward -- an access violation is 0xC0000005 =
    3221225477. That is why 1.0's abort on a locus with no observed allele
    reads as a signal on Linux and as a large positive status on Windows.

    The top 256 values are excluded deliberately: 4294967295 is exit(-1) seen
    as unsigned, an ordinary failure, and reading it as a crash made three
    one_group cases -- where both builds correctly refuse k > 1 for a single
    grouping -- fail on Windows while passing everywhere else.
    """
    return returncode < 0 or 0xC0000000 <= returncode < 0xFFFFFF00


def is_results_file(name):
    return not (name.endswith("_summary") or name.endswith(".summary.txt")
                or name.endswith("deletedloci"))


def said_about_itself(directory, width=120):
    """The last thing a run printed, for a file it wrote nothing into.

    An empty result file is not a numerical disagreement, it is a run that
    declined to report -- and in the legacy layout an undefined row is omitted
    rather than marked, so the file goes empty without anything in it saying
    why. What the run printed is where the reason is.
    """
    for name in ("console.err", "console.out"):
        path = os.path.join(directory, name)
        if not os.path.exists(path):
            continue
        said = [l.strip() for l in open(path, errors="replace") if l.strip()]
        if said:
            return "; it said: %r" % said[-1][:width]
    return ""


def first_difference(ref_bytes, cnd_bytes, width=64):
    """Where two masked files first disagree, short enough to sit in a log.

    A bare list of differing filenames says nothing about the cause -- on one
    machine it is a number, on another it could be a header, a row that should
    not be there, or a line ending. Printing the first disagreeing line makes
    the log self-diagnosing, which matters most when the failure is on a
    machine you cannot reach.
    """
    ref = ref_bytes.decode("utf-8", "replace").splitlines()
    cnd = cnd_bytes.decode("utf-8", "replace").splitlines()

    for i in range(max(len(ref), len(cnd))):
        r = ref[i] if i < len(ref) else None
        c = cnd[i] if i < len(cnd) else None
        if r == c:
            continue
        if r is None:
            return "line %d: candidate has %r, reference ends" % (i + 1, c[:width])
        if c is None:
            return "line %d: reference has %r, candidate ends" % (i + 1, r[:width])
        return "line %d: ref %r vs cnd %r" % (i + 1, r[:width], c[:width])

    return "%d vs %d bytes, no differing line" % (len(ref_bytes), len(cnd_bytes))


def ref_fulldata_cells(path):
    """{(label, g, locus): printed value} from version 1.0's _fulldata layout.

    1.0 wrote one row per label and g, with a column per locus: the locus
    names are in the header, between the label columns and MEAN VAR STD_ERR.
    Tuple files carry one label column per tuple member; they are joined with
    a comma to match the single label column the candidate writes.
    """
    lines = open(path).read().splitlines()
    if not lines:
        return {}
    head = lines[0].split()
    if "G" not in head or "NUM_LOCI" not in head:
        return {}
    k = head.index("G")                     # label columns before it
    names = head[k + 2:-3]                  # ... G NUM_LOCI <loci> MEAN VAR STD_ERR
    cells = {}
    for line in lines[1:]:
        f = line.split()
        if len(f) < k + 3 + len(names):
            continue
        label = ",".join(f[:k])
        g = f[k]
        for name, value in zip(names, f[k + 2:k + 2 + len(names)]):
            cells[(label, g, name)] = value
    return cells


def cnd_fulldata_cells(path):
    """The same, from the candidate's layout: one row per label and locus."""
    lines = open(path).read().splitlines()
    if not lines:
        return {}
    head = lines[0].split("\t")
    if len(head) < 3 or head[1] != "LOCUS":
        return {}
    gs = [h[1:] for h in head[2:]]          # G1, G2, ... -> 1, 2, ...
    cells = {}
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) != len(head):
            continue
        for g, value in zip(gs, f[2:]):
            cells[(f[0], g, f[1])] = value
    return cells


def compare_fulldata(a, b):
    """Every per-locus value 1.0 printed, found unchanged in the candidate.

    The two layouts are transposes of each other (see the manual, "Changes
    from version 1.0"), so this compares cells rather than bytes. The
    candidate holds cells the reference does not -- g = 1, and the rows 1.0
    dropped because some locus was undefined -- which are not evidence of
    disagreement and are not compared here.
    """
    ref, cnd = ref_fulldata_cells(a), cnd_fulldata_cells(b)
    if not ref:
        return "reference _fulldata has no rows to compare"
    absent = [k for k in ref if k not in cnd]
    if absent:
        k = sorted(absent)[0]
        return "%d of %d cells absent from the candidate, e.g. %s g=%s %s" % (
            len(absent), len(ref), k[0], k[1], k[2])
    bad = [k for k in ref if ref[k] != cnd[k]]
    if bad:
        k = sorted(bad)[0]
        return "%d of %d cells differ, e.g. %s g=%s %s: %s vs %s" % (
            len(bad), len(ref), k[0], k[1], k[2], ref[k], cnd[k])
    return None


#Summary columns: how close is "the same number".
#
#The candidate forms the mean and variance in one pass (Welford) where 1.0
#summed twice over the stored per-locus values, so the last bits of a variance
#can differ. SUMMARY_RTOL is well below the six significant digits that get
#printed, so a genuine disagreement still fails; SUMMARY_ATOL is what lets
#"both are zero to within rounding" pass -- 1.0 printed 2.96e-31 for a
#variance whose every deviation was zero, where the one-pass form gives 0.
SUMMARY_RTOL = 1e-9
SUMMARY_ATOL = 1e-12
SUMMARY_COLUMNS = 3


def same_number(x, y):
    if x == y:
        return True
    try:
        a, b = float(x), float(y)
    except ValueError:
        return False
    if a != a and b != b:                      # both nan
        return True
    if abs(a) <= SUMMARY_ATOL and abs(b) <= SUMMARY_ATOL:
        return True
    return abs(a - b) <= SUMMARY_RTOL * max(abs(a), abs(b))


def compare_stats(ref_bytes, cnd_bytes):
    """Rows of a statistics file: labels exact, the summary to tolerance."""
    ref = [l for l in ref_bytes.decode().splitlines() if l.strip()]
    cnd = [l for l in cnd_bytes.decode().splitlines() if l.strip()]
    if len(ref) != len(cnd):
        return "%d rows vs %d" % (len(ref), len(cnd))

    for i, (a, b) in enumerate(zip(ref, cnd)):
        fa, fb = a.split(), b.split()
        if len(fa) != len(fb) or fa[:-SUMMARY_COLUMNS] != fb[:-SUMMARY_COLUMNS]:
            return "line %d: ref %r vs cnd %r" % (i + 1, a[:64], b[:64])
        for x, y in zip(fa[-SUMMARY_COLUMNS:], fb[-SUMMARY_COLUMNS:]):
            if not same_number(x, y):
                return "line %d: ref %r vs cnd %r" % (i + 1, a[:64], b[:64])
    return None


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

        if f.endswith("_fulldata"):
            why = compare_fulldata(a, b)
            if why:
                diffs.append("%s (%s)" % (f, why))
            continue

        ref_bytes = masked(a)
        cnd_bytes = masked(b, drop_g1=is_results_file(f))

        if is_results_file(f) and ref_bytes.strip() and cnd_bytes.strip():
            why = compare_stats(ref_bytes, cnd_bytes)
            if why:
                diffs.append("%s (%s)" % (f, why))
            else:
                n, labels = g1_rows(b)
                if n != len(labels):
                    missing.append("%s: %d g=1 rows for %d labels" % (f, n, len(labels)))
            continue

        if ref_bytes != cnd_bytes:
            if not ref_bytes.strip() and cnd_bytes.strip():
                why = "reference wrote no rows%s" % said_about_itself(ref_dir)
            elif not cnd_bytes.strip() and ref_bytes.strip():
                why = "candidate wrote no rows%s" % said_about_itself(cnd_dir)
            else:
                why = first_difference(ref_bytes, cnd_bytes)
            diffs.append("%s (%s)" % (f, why))
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
        # 2.0 writes tab-separated output with a header by default, so the
        # candidate is asked for 1.0's layout explicitly. The reference cannot
        # be given the flag -- it predates it -- which is why the two argv
        # lists differ here and nowhere else.
        cnd = run_one(args.candidate, cnd_dir, dfile, meta, "out", kw,
                      extra_argv=("--legacy",))

        if died(ref.returncode):
            # The reference died: nothing to compare against.
            if died(cnd.returncode):
                status, detail = "XFAIL", (
                    "reference and candidate both die (status %d) "
                    "-- pre-existing defect" % ref.returncode)
            elif cnd.returncode != 0:
                status, detail = "FAIL", (
                    "reference died (status %d), candidate exited %d" % (
                        ref.returncode, cnd.returncode))
            else:
                status, detail = "FIXED", (
                    "reference died (status %d), candidate exited 0 with %d outputs" % (
                        ref.returncode, len(outputs(cnd_dir, dfile))))
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
