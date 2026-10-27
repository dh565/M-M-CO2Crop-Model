#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial
# research use with citation; redistribution of modified versions and commercial use
# require the copyright holders' consent.
"""BLOCKING release-acceptance test for the M&M-CO2Wheat.v1.0 code clean-up.

The test it encodes, in one sentence: *the current tree compiles and reproduces the sealed
v1.0 reference run package byte for byte on all 42 decks.*  Nothing else is an acceptance criterion --
not a scorecard, not a plot, not a prose claim.  A clean-up commit that changes one
output byte on one deck fails here and is not a clean-up.

Why the REFERENCE ARM's OWN `input/` folders are the inputs, never the tree's
`input/decks`: the reference package was assembled by a deck assembler that BUILDS each
deck's `driver_main.txt` (the deck's own on-disk lines + the reference template's tail)
and resolves some manifest rows to the then-current root `input/`.  Re-running that
assembly would re-derive the inputs from a moving tree and so would test the assembler,
not the model.  The sealed arm already carries the exact bytes the frozen executable was
fed, deck by deck, so this harness copies those bytes verbatim and asks the only question
a code clean-up can be asked: given identical inputs, does the rebuilt executable write
identical outputs?

Invocation convention, taken from the arm itself (scripts/run_decks_sites.py:run_one and
the arm's deck folders, which hold input/, output/ and run.log and NO exe copy):
  - the executable is invoked by ABSOLUTE path, `cwd` = the deck folder;
  - `output/` is created before the run (the model does not create it);
  - stdout and stderr both go to `run.log`, and the runner then appends `exit=<rc>\n`
    to the SAME open handle, with no separating newline (a one-byte artifact);
  - liveness is checked on the child's own PID, never a name-based pgrep sweep (this host
    runs other sessions' copies of the same binary name concurrently).

Comparison:
  - every file under <out>/<deck>/output vs the arm's <deck>/output, byte comparison
    (filecmp.cmp shallow=False), classified identical / differing / missing / extra;
  - run.log compared separately and reported separately (it carries the banner text and
    the season summaries, which are model output, but is not part of the identity verdict
    because a runner is free to add wall-clock lines; a raw-byte match is reported when
    it holds, and a timestamp-stripped comparison is the fallback);
  - the copied input/ is md5-verified against the arm's, so an input-copy failure can
    never be mistaken for a model difference.

Verdict: PASS only if all 42 decks exit 0, no deck left a live process, every deck's
input copy verified, and the totals are 0 differing, 0 missing and 0 extra output files.

Usage:
  release_acceptance.py --exe <mmco2crop> --ref-arm <reference package dir> \
      --out <clock-named output dir> [--workers N] [--decks A,B,...] [--label TEXT]
"""
import argparse
import concurrent.futures
import filecmp
import hashlib
import os
import shutil
import subprocess
import sys
import time

TS_HINTS = ("Mon ", "Tue ", "Wed ", "Thu ", "Fri ", "Sat ", "Sun ",
            "elapsed", "Elapsed", "wall", "Wall", "sec)", "seconds")


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel_files(root):
    """Sorted relative paths of every regular file under root ('' if root absent)."""
    out = []
    if not os.path.isdir(root):
        return out
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            p = os.path.join(dirpath, f)
            if os.path.isfile(p):
                out.append(os.path.relpath(p, root))
    return sorted(out)


def pid_alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def first_diff_line(a, b):
    """1-indexed first differing line of two text files, with both sides truncated."""
    try:
        with open(a, "r", errors="replace") as fa, open(b, "r", errors="replace") as fb:
            i = 0
            while True:
                la, lb = fa.readline(), fb.readline()
                i += 1
                if not la and not lb:
                    return None
                if la != lb:
                    return (i, la.rstrip("\n")[:300], lb.rstrip("\n")[:300])
    except (OSError, UnicodeDecodeError):
        return (0, "<binary or unreadable>", "<binary or unreadable>")


def strip_ts(path):
    """run.log lines with no wall-clock hint -- the fallback run.log comparison."""
    keep = []
    with open(path, "r", errors="replace") as fh:
        for ln in fh:
            if any(h in ln for h in TS_HINTS):
                continue
            keep.append(ln)
    return keep


def run_one(deck, ref_arm, out, exe):
    """Stage the arm's own input/ bytes, run the exe deck-local, verify the input copy."""
    ref_deck = os.path.join(ref_arm, deck)
    ddir = os.path.join(out, deck)
    if os.path.exists(ddir):
        shutil.rmtree(ddir)
    os.makedirs(ddir)
    shutil.copytree(os.path.join(ref_deck, "input"), os.path.join(ddir, "input"))
    # Input-copy verification: a staging failure must never be read as a model difference.
    bad_in = []
    for rp in rel_files(os.path.join(ref_deck, "input")):
        src = os.path.join(ref_deck, "input", rp)
        dst = os.path.join(ddir, "input", rp)
        if not os.path.isfile(dst) or md5_file(src) != md5_file(dst):
            bad_in.append(rp)
    os.makedirs(os.path.join(ddir, "output"), exist_ok=True)
    t0 = time.time()
    log_path = os.path.join(ddir, "run.log")
    with open(log_path, "w") as logf:
        proc = subprocess.Popen([os.path.abspath(exe)], cwd=ddir,
                                stdout=logf, stderr=subprocess.STDOUT)
        rc = proc.wait()
        logf.write("exit=%d\n" % rc)
    sec = time.time() - t0
    time.sleep(0.2)
    live = 1 if pid_alive(proc.pid) else 0
    return {"deck": deck, "sec": sec, "exit": rc, "live": live, "bad_input": bad_in}


def compare_deck(deck, ref_arm, out):
    ref_out = os.path.join(ref_arm, deck, "output")
    new_out = os.path.join(out, deck, "output")
    ref_set, new_set = rel_files(ref_out), rel_files(new_out)
    missing = sorted(set(ref_set) - set(new_set))
    extra = sorted(set(new_set) - set(ref_set))
    same, diff, diff_detail = [], [], []
    for rp in sorted(set(ref_set) & set(new_set)):
        a, b = os.path.join(ref_out, rp), os.path.join(new_out, rp)
        if filecmp.cmp(a, b, shallow=False):
            same.append(rp)
        else:
            diff.append(rp)
            if len(diff_detail) < 5:
                diff_detail.append((rp, first_diff_line(a, b)))
    # run.log, reported but NOT part of the identity verdict
    rl_ref = os.path.join(ref_arm, deck, "run.log")
    rl_new = os.path.join(out, deck, "run.log")
    if os.path.isfile(rl_ref) and os.path.isfile(rl_new):
        if filecmp.cmp(rl_ref, rl_new, shallow=False):
            rl = "IDENTICAL"
        elif strip_ts(rl_ref) == strip_ts(rl_new):
            rl = "IDENTICAL_AFTER_TS_STRIP"
        else:
            fd = first_diff_line(rl_ref, rl_new)
            rl = "DIFFERS@line%s" % (fd[0] if fd else "?")
    else:
        rl = "MISSING"
    return {"deck": deck, "identical": len(same), "differing": len(diff),
            "missing": len(missing), "extra": len(extra), "runlog": rl,
            "diff_files": diff, "missing_files": missing, "extra_files": extra,
            "diff_detail": diff_detail}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True)
    ap.add_argument("--ref-arm", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--decks", default="", help="comma-separated subset (default: all in the arm)")
    ap.add_argument("--label", default="", help="free-text line recorded in ACCEPTANCE_LOG.md")
    ap.add_argument("--flags", default="", help="compiler flags, recorded verbatim in the log")
    args = ap.parse_args()

    exe = os.path.abspath(args.exe)
    ref_arm = os.path.abspath(args.ref_arm)
    out = os.path.abspath(args.out)
    if not os.path.isfile(exe):
        sys.exit("ERROR: executable not found: %s" % exe)
    if not os.path.isdir(ref_arm):
        sys.exit("ERROR: reference arm not found: %s" % ref_arm)
    decks = sorted(d for d in os.listdir(ref_arm)
                   if os.path.isdir(os.path.join(ref_arm, d, "input"))
                   and os.path.isdir(os.path.join(ref_arm, d, "output")))
    if args.decks:
        want = args.decks.split(",")
        missing = [d for d in want if d not in decks]
        if missing:
            sys.exit("ERROR: decks not in arm: %s" % missing)
        decks = want
    os.makedirs(out, exist_ok=True)
    exe_md5 = md5_file(exe)
    t_start = time.time()
    print("exe       %s" % exe)
    print("exe md5   %s" % exe_md5)
    print("ref arm   %s  (%d decks)" % (ref_arm, len(decks)))
    print("out       %s" % out)
    print("%-24s %8s %6s %6s" % ("deck", "sec", "exit", "live"))
    runs = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(run_one, d, ref_arm, out, exe): d for d in decks}
        for fut in concurrent.futures.as_completed(futs):
            r = fut.result()
            runs[r["deck"]] = r
            print("%-24s %8.1f %6d %6d%s" % (r["deck"], r["sec"], r["exit"], r["live"],
                  "  BAD_INPUT_COPY=%d" % len(r["bad_input"]) if r["bad_input"] else ""))
            sys.stdout.flush()

    cmps = {d: compare_deck(d, ref_arm, out) for d in decks}
    wall = time.time() - t_start

    tot = {k: sum(cmps[d][k] for d in decks) for k in ("identical", "differing", "missing", "extra")}
    n_exit0 = sum(1 for d in decks if runs[d]["exit"] == 0)
    n_live = sum(runs[d]["live"] for d in decks)
    n_badin = sum(len(runs[d]["bad_input"]) for d in decks)
    rl_ident = sum(1 for d in decks if cmps[d]["runlog"] == "IDENTICAL")
    rl_strip = sum(1 for d in decks if cmps[d]["runlog"] == "IDENTICAL_AFTER_TS_STRIP")
    expect_42 = len(decks) == 42
    verdict = ("PASS" if (n_exit0 == len(decks) and n_live == 0 and n_badin == 0
                          and tot["differing"] == 0 and tot["missing"] == 0
                          and tot["extra"] == 0 and expect_42) else "FAIL")

    # MD5SUMS of everything this run produced (run.log + output/), md5sum-compatible.
    md5_path = os.path.join(out, "MD5SUMS")
    n_md5 = 0
    with open(md5_path, "w") as fh:
        for d in decks:
            ddir = os.path.join(out, d)
            for rp in ["run.log"] + [os.path.join("output", x)
                                     for x in rel_files(os.path.join(ddir, "output"))]:
                p = os.path.join(ddir, rp)
                if os.path.isfile(p):
                    fh.write("%s  %s/%s\n" % (md5_file(p), d, rp))
                    n_md5 += 1

    stamp = time.strftime("%Y-%m-%d %H:%M:%S %Z", time.localtime())
    L = []
    L.append("# Release acceptance -- v1.0 code clean-up identity test")
    L.append("")
    L.append("**The test:** the current tree compiles and reproduces the sealed v1.0 reference run package")
    L.append("byte for byte on all %d decks. PASS is the only acceptable result for a" % len(decks))
    L.append("clean-up change; a single differing output byte on a single deck is a FAIL.")
    L.append("")
    L.append("| field | value |")
    L.append("|---|---|")
    L.append("| verdict | **%s** |" % verdict)
    L.append("| date | %s |" % stamp)
    L.append("| label | %s |" % (args.label or "(none)"))
    L.append("| executable | `%s` |" % exe)
    L.append("| exe md5 | `%s` |" % exe_md5)
    L.append("| compiler flags | `%s` |" % (args.flags or "(not recorded)"))
    L.append("| reference arm | `%s` |" % ref_arm)
    L.append("| output | `%s` |" % out)
    L.append("| decks | %d (expected 42: %s) |" % (len(decks), "yes" if expect_42 else "NO"))
    L.append("| workers | %d |" % args.workers)
    L.append("| wall time (run + compare) | %.1f s |" % wall)
    L.append("| decks exit 0 | %d / %d |" % (n_exit0, len(decks)))
    L.append("| live processes left | %d |" % n_live)
    L.append("| input files failing md5 copy check | %d |" % n_badin)
    L.append("| output files identical | %d |" % tot["identical"])
    L.append("| output files differing | %d |" % tot["differing"])
    L.append("| output files missing | %d |" % tot["missing"])
    L.append("| output files extra | %d |" % tot["extra"])
    L.append("| run.log identical (raw bytes) | %d / %d |" % (rl_ident, len(decks)))
    L.append("| run.log identical after timestamp strip only | %d / %d |" % (rl_strip, len(decks)))
    L.append("| MD5SUMS entries written | %d |" % n_md5)
    L.append("")
    L.append("run.log is reported separately and is NOT part of the identity verdict: a runner")
    L.append("may add wall-clock lines a model-identity test must not fail on. Its banner and")
    L.append("season-summary content is model output, so a raw-byte run.log match is the")
    L.append("stronger signal and is recorded above when it holds.")
    L.append("")
    L.append("## Per-deck table")
    L.append("")
    L.append("| deck | sec | exit | live | identical | differing | missing | extra | run.log |")
    L.append("|---|---:|---:|---:|---:|---:|---:|---:|---|")
    for d in decks:
        r, c = runs[d], cmps[d]
        L.append("| %s | %.1f | %d | %d | %d | %d | %d | %d | %s |"
                 % (d, r["sec"], r["exit"], r["live"], c["identical"], c["differing"],
                    c["missing"], c["extra"], c["runlog"]))
    L.append("")
    if tot["differing"] or tot["missing"] or tot["extra"] or n_badin:
        L.append("## Non-identity detail")
        L.append("")
        for d in decks:
            c = cmps[d]
            if not (c["differing"] or c["missing"] or c["extra"] or runs[d]["bad_input"]):
                continue
            L.append("### %s" % d)
            if runs[d]["bad_input"]:
                L.append("- input copy failed md5: %s" % ", ".join(runs[d]["bad_input"][:20]))
            if c["missing_files"]:
                L.append("- missing (%d): %s" % (len(c["missing_files"]), ", ".join(c["missing_files"][:20])))
            if c["extra_files"]:
                L.append("- extra (%d): %s" % (len(c["extra_files"]), ", ".join(c["extra_files"][:20])))
            if c["diff_files"]:
                L.append("- differing (%d): %s" % (len(c["diff_files"]), ", ".join(c["diff_files"][:20])))
                for rp, fd in c["diff_detail"]:
                    if fd:
                        L.append("  - `%s` first differs at line %d" % (rp, fd[0]))
                        L.append("    - reference: `%s`" % fd[1])
                        L.append("    - new: `%s`" % fd[2])
            L.append("")
    L.append("Produced by `scripts/release_acceptance.py`.")
    with open(os.path.join(out, "ACCEPTANCE_LOG.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")

    print("")
    print("VERDICT=%s  exit0=%d/%d live=%d bad_input=%d identical=%d differing=%d missing=%d extra=%d"
          % (verdict, n_exit0, len(decks), n_live, n_badin,
             tot["identical"], tot["differing"], tot["missing"], tot["extra"]))
    print("run.log raw-identical %d/%d, ts-strip-identical %d/%d" % (rl_ident, len(decks), rl_strip, len(decks)))
    print("log  %s" % os.path.join(out, "ACCEPTANCE_LOG.md"))
    print("md5  %s (%d entries)" % (md5_path, n_md5))
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
