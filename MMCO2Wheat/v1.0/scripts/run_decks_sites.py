#!/usr/bin/env python3
"""Three-site BASELINE runner: harmonizes Braunschweig / Maricopa / AgMIP / WDF23 / WDF25
site decks onto the CURRENT, FULL v1.0 driver template and runs them.

Generalizes `run_decks_rec130.py` (rewrite-in-place discipline, md5/exit=/live-process
checks) to a DIFFERENT problem: deck families carry `driver_main.txt` files of DIFFERENT,
shorter lengths than the current reference template -- on disk today: Maricopa/
Braunschweig 116 lines, AgMIP 144 lines, WDF23 118 lines -- each built before some later
record range existed. Every record beyond a deck's own on-disk length is a crop-MECHANISM
switch (not a site input: see the reference driver's own comments at those lines), so "one
executable, one stack, all sites" means APPENDING the template's own tail verbatim from the
reference driver onto every site deck's own (unmodified) file, never renumbering or
re-deriving anything, and NEVER hard-coding one specific deck length or template length.

2026-10-11 generalization (the records 145/146 adoption, records 145/146): the previous version
of this script hard-coded deck length == 116 and reference length == 141 (`ValueError`
otherwise), so it could not stage a deck past record 141 at all, and could not reach a
146-line template (records 145/146) on ANY deck. `build_driver_main()` below now reads
BOTH lengths from the files themselves and appends `ref_lines[len(deck_lines):]` -- the
template's own missing tail, whatever its length -- so a 116-, 118- or 144-line on-disk
deck driver is staged against a 141-, 144- or 146-line (or any future length) template with
no code change, and the deck's own on-disk lines are copied verbatim, never edited (so they
stay byte-exact on every overlapping line) except for an explicit, logged `--force-rec`.

Per-deck assembly (manifest-driven, generalizing `deck_driver_gen.py`'s ROOT|DECK
resolution to the Braunschweig manifest as well as the Maricopa/AgMIP/WDF25 one):
  for each (deck, file) row in the given DECK_CONFIG_MANIFEST.csv:
    source = <root_input_dir>/<file>   if row['source'] == 'ROOT'
             <decks_dir>/<deck>/<file> otherwise (== 'DECK')
  copied verbatim into <out>/<deck>/input/<file> -- EXCEPT driver_main.txt, which is
  built as: <decks_dir>/<deck>/driver_main.txt lines 1-N (the deck's own on-disk length N,
  verbatim, byte-for-byte, never asserted to any one fixed N) + <ref_driver> lines N+1 to
  the template's own length (verbatim, including every comment). This also means
  driver_minerals.csv and driver_resp_mainten.txt resolve to the CURRENT root input/
  copies (ROOT in both manifests), not any stale deck-local snapshot -- the harmonization
  log records old-vs-new md5/column-count for every such substitution so a stale
  15-numeric-column driver_minerals.csv (missing the rec-131 avail_A_mgL_CAL /
  xylem_Ca_mgL_CAL columns) is replaced by the current 18-column file and the substitution
  is machine-visible, not asserted in prose.

Usage:
  run_decks_sites.py <manifest.csv> <decks_dir> <root_input_dir> <ref_driver> <out_pkg> \
      --exe <path to mmco2crop> --decks NAME,NAME,... [--workers N]

Exit 0 only if every requested deck exists in the manifest, every copy/append
succeeded, every run returned exit code 0, and no live mmco2crop process was
observed after any run returned.
"""
import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d6_config_guard  # noqa: E402  -- D6 2026-10-10 configuration-harmonization guard


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_manifest(path):
    decks = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            decks.setdefault(row["deck"], []).append(row)
    return decks


def trailing_numeric_cols(path):
    """Column count of the last non-'!'-comment, non-blank line, comma-split."""
    last = None
    with open(path) as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith("!"):
                continue
            last = s
    if last is None:
        return 0
    return len(last.split(","))


def build_driver_main(deck_driver_path, ref_driver_path, out_path, log_rows, deck, force_recs=None):
    """force_recs: {1-indexed line number: value} -- a DECLARED, LOGGED override of a
    record inside the site's OWN on-disk range, used ONLY for a record this task found to
    be a hard ERROR STOP blocker when the template's own tail is appended as-is (e.g.
    record 99 tiller_survival_mech, required == 4 by the appended record 123
    tiller_n_afford=1). This is the SAME 'minimal documented fix, marked
    DECLARED-TRANSFERRED' the brief authorizes for a missing auxiliary input, applied to a
    record instead of a file, because the record is itself a crop-MECHANISM switch (its
    own in-source comment names it as such), never a site measurement -- and is reported,
    never silent.

    2026-10-11 generalization (the records 145/146 adoption): NEITHER the deck's own length NOR
    the reference template's length is hard-coded (the previous version asserted exactly
    116 and exactly 141, which made it impossible to stage a deck past record 141, or to
    reach a 146-line template on any deck, without editing this function). Both lengths
    are read from the files themselves; the template's own tail, `ref_lines[n_deck:]`,
    whatever its length, supplies every record the deck does not itself carry. The deck's
    own on-disk lines are never rewritten (so they stay byte-exact on every overlapping
    line) except for an explicit, logged `--force-rec` on a line within the deck's own
    length -- a `--force-rec` targeting a line beyond the deck's own length is refused
    (there is nothing there to override; that record only exists in the appended tail,
    which is already the template's own value by construction)."""
    force_recs = force_recs or {}
    with open(deck_driver_path) as fh:
        deck_lines = [ln.rstrip("\n") for ln in fh]
    with open(ref_driver_path) as fh:
        ref_lines = [ln.rstrip("\n") for ln in fh]
    n_deck = len(deck_lines)
    n_ref = len(ref_lines)
    if n_deck > n_ref:
        raise ValueError(
            "deck %s driver_main.txt has %d lines, longer than the reference template's %d"
            " lines -- the template must supply every line beyond the deck's own; a deck"
            " longer than the template cannot be staged by appending a tail"
            % (deck, n_deck, n_ref)
        )
    for lineno, value in sorted(force_recs.items()):
        if lineno < 1 or lineno > n_deck:
            raise ValueError(
                "deck %s: --force-rec %d:%d targets a line outside the deck's own 1-%d"
                " on-disk range (force-rec only overrides a line the deck itself carries;"
                " a line beyond %d is already the template's own appended value)"
                % (deck, lineno, value, n_deck, n_deck)
            )
        idx = lineno - 1
        old_tok = deck_lines[idx].split()[0] if deck_lines[idx].split() else "<blank>"
        deck_lines[idx] = "%d    ! rec%d FORCED to the reference stack's value by run_decks_sites.py" \
            " (DECLARED-TRANSFERRED mechanism switch, site built value was %s; see HARMONIZATION_LOG.csv)" \
            % (value, lineno, old_tok)
        log_rows.append((deck, "driver_main.txt", "line_%d" % lineno,
                          "FORCED_SITE_RECORD old=%s new=%s (blocking ERROR STOP dependency)" % (old_tok, value), str(value)))
    appended = ref_lines[n_deck:n_ref]  # the template's own tail, whatever the deck's own length
    out_lines = deck_lines + appended
    with open(out_path, "w") as fh:
        fh.write("\n".join(out_lines) + "\n")
    # Independent verification: re-read the WRITTEN file (not the in-memory lists) and
    # confirm line count, that every ORIGINAL deck line is byte-exact (except a declared
    # --force-rec line, which is expected to differ from the on-disk deck file but whose
    # leading token must equal the forced value), and that every appended record's leading
    # token matches the reference's, line by line.
    with open(out_path) as fh:
        written = [ln.rstrip("\n") for ln in fh]
    if len(written) != n_ref:
        raise ValueError("deck %s: written driver_main.txt has %d lines, expected %d (template length)"
                          % (deck, len(written), n_ref))
    with open(deck_driver_path) as fh:
        deck_lines_ondisk = [ln.rstrip("\n") for ln in fh]
    for i in range(n_deck):
        lineno = i + 1
        if lineno in force_recs:
            got_tok = written[i].split()[0]
            want_tok = str(force_recs[lineno])
            if got_tok != want_tok:
                raise ValueError(
                    "deck %s: forced line %d token mismatch after write: got %r want %r"
                    % (deck, lineno, got_tok, want_tok)
                )
        elif written[i] != deck_lines_ondisk[i]:
            raise ValueError(
                "deck %s: line %d (within the deck's own on-disk range) is not byte-exact"
                " after staging -- the deck's own lines must be copied verbatim unless"
                " declared via --force-rec" % (deck, lineno)
            )
    for i in range(n_deck, n_ref):
        got_tok = written[i].split()[0]
        want_tok = ref_lines[i].split()[0]
        if got_tok != want_tok:
            raise ValueError(
                "deck %s: line %d token mismatch after append: got %r want %r"
                % (deck, i + 1, got_tok, want_tok)
            )
        log_rows.append((deck, "driver_main.txt", "line_%d" % (i + 1), "APPENDED_FROM_REF", got_tok))
    return n_deck, len(appended)


def assemble_deck(deck, rows, decks_dir, root_input, ref_driver, out_deck_dir, log_rows, force_recs=None):
    in_dir = os.path.join(out_deck_dir, "input")
    os.makedirs(in_dir, exist_ok=True)
    for row in rows:
        fname = row["file"]
        source = row["source"]
        src_path = (os.path.join(root_input, fname) if source == "ROOT"
                    else os.path.join(decks_dir, deck, fname))
        dst_path = os.path.join(in_dir, fname)
        if fname == "driver_main.txt":
            continue  # handled specially below, after the plain copies
        if not os.path.isfile(src_path):
            raise ValueError("deck %s: manifest file %s (%s) not found at %s" % (deck, fname, source, src_path))
        # D6 2026-10-10 bug fix: a manifest row may name a file under a subdirectory
        # (e.g. "WDF25_Data/foo.csv" -- the first manifest this script was run against
        # with such a path); the previous code only created the top-level `input/`
        # directory, so the destination copy failed with FileNotFoundError on any
        # nested relative path. Create the destination's own parent directory first.
        os.makedirs(os.path.dirname(dst_path), exist_ok=True)
        shutil.copyfile(src_path, dst_path)
        old_md5 = row.get("md5", "")
        new_md5 = md5_file(dst_path)
        changed = "CHANGED" if new_md5 != old_md5 else "unchanged"
        extra = ""
        if fname in ("driver_minerals.csv",):
            extra = "ncols=%d" % trailing_numeric_cols(dst_path)
        log_rows.append((deck, fname, source, "manifest_md5=%s run_md5=%s %s (%s)"
                          % (old_md5, new_md5, extra, changed), ""))
    # driver_main.txt: site lines 1-N (the deck's OWN on-disk length N, from decks_dir, DECK by
    # construction, verbatim) + the reference template's own tail, lines N+1 to the template's
    # own length (verbatim) -- N and the template length are read from the files, never assumed.
    deck_driver_rows = [r for r in rows if r["file"] == "driver_main.txt"]
    if not deck_driver_rows:
        raise ValueError("deck %s: manifest carries no driver_main.txt row" % deck)
    deck_driver_src = os.path.join(decks_dir, deck, "driver_main.txt")
    out_driver = os.path.join(in_dir, "driver_main.txt")
    n_site, n_appended = build_driver_main(deck_driver_src, ref_driver, out_driver, log_rows, deck, force_recs)
    log_rows.append((deck, "driver_main.txt", "ASSEMBLY",
                      "site_lines=%d (1-%d, verbatim) + ref_lines=%d (%d-%d, verbatim)"
                      % (n_site, n_site, n_appended, n_site + 1, n_site + n_appended), ""))


def pid_alive(pid):
    """Exact-PID liveness check (never a name-based pgrep sweep: this host runs OTHER
    sessions' copies of the SAME binary name concurrently -- confirmed live during this
    task's own smoke test -- so `pgrep -x mmco2crop` reports processes that are
    not this run's child and is not a valid 'did MY run leave something behind' test.
    `subprocess.run` already waitpid()s its own child to exit before returning, so the
    only question worth asking is whether THIS PID, specifically, still exists)."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # exists, owned by someone else -- should not happen for our own child
    return True


def run_one(deck, out_deck_dir, exe, exe_basename):
    os.makedirs(os.path.join(out_deck_dir, "output"), exist_ok=True)
    t0 = time.time()
    log_path = os.path.join(out_deck_dir, "run.log")
    with open(log_path, "w") as logf:
        proc = subprocess.Popen([exe], cwd=out_deck_dir, stdout=logf, stderr=subprocess.STDOUT)
        rc = proc.wait()
        logf.write("exit=%d\n" % rc)
    dt = time.time() - t0
    time.sleep(0.2)
    n_live = 1 if pid_alive(proc.pid) else 0
    return deck, dt, rc, n_live


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest")
    ap.add_argument("decks_dir")
    ap.add_argument("root_input")
    ap.add_argument("ref_driver")
    ap.add_argument("out_pkg")
    ap.add_argument("--exe", required=True)
    ap.add_argument("--decks", required=True, help="comma-separated deck subset")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--force-rec", action="append", default=[],
                     help="LINE:VALUE (1-indexed driver_main.txt line, within the deck's OWN"
                          " on-disk length -- refused otherwise) to force to the reference"
                          " stack's value -- a DECLARED, LOGGED override for a record this"
                          " task found to hard-ERROR-STOP under the appended template tail as"
                          " built for the site (e.g. --force-rec 99:4). Repeatable.")
    ap.add_argument("--config-reference", default=d6_config_guard.DEFAULT_REFERENCE,
                     help="D6 2026-10-10 guard: path to the canonical v1.0 driver_main.txt"
                          " (default docs/V1_CONFIGURATION_driver_main_v1.0.txt) every assembled"
                          " deck's MECHANISM records (23 to the reference file's own length,"
                          " excluding the declared SITE-STATE/DIAGNOSTIC records, see"
                          " scripts/d6_config_guard.py) are checked against before any run is"
                          " launched.")
    ap.add_argument("--allow-config-deviation", action="store_true",
                     help="D6 2026-10-10 guard: do not refuse on a MECHANISM-record mismatch against"
                          " --config-reference; log every deviation to CONFIG_DEVIATION_LOG.csv in"
                          " the output package instead and proceed. For a DECLARED ablation arm only"
                          " -- never for a routine site/fingerprint run.")
    args = ap.parse_args()
    force_recs = {}
    for spec in args.force_rec:
        lineno_s, value_s = spec.split(":")
        force_recs[int(lineno_s)] = int(value_s)

    manifest = os.path.abspath(args.manifest)
    decks_dir = os.path.abspath(args.decks_dir)
    root_input = os.path.abspath(args.root_input)
    ref_driver = os.path.abspath(args.ref_driver)
    out = os.path.abspath(args.out_pkg)
    exe = os.path.abspath(args.exe)
    exe_basename = os.path.basename(exe)
    if not os.path.isfile(exe):
        print("ERROR: executable not found: %s" % exe); sys.exit(2)
    os.makedirs(out, exist_ok=True)

    manifest_decks = load_manifest(manifest)
    wanted = args.decks.split(",")
    missing = [d for d in wanted if d not in manifest_decks]
    if missing:
        print("ERROR: decks not in manifest %s: %s" % (manifest, missing)); sys.exit(2)

    log_rows = []
    out_dirs = {}
    for d in wanted:
        ddir = os.path.join(out, d)
        if os.path.exists(ddir):
            shutil.rmtree(ddir)
        assemble_deck(d, manifest_decks[d], decks_dir, root_input, ref_driver, ddir, log_rows, force_recs)
        out_dirs[d] = ddir

    with open(os.path.join(out, "HARMONIZATION_LOG.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["deck", "file", "source_or_marker", "detail", "token"])
        w.writerows(log_rows)

    # D6 2026-10-10 configuration-harmonization guard: every assembled deck's MECHANISM
    # records must match the canonical v1.0 reference before any run is launched. A
    # record FORCED by --force-rec is, by construction, declared and logged above (not a
    # silent deviation) and is skipped here only if it was actually forced this run.
    config_reference = os.path.abspath(args.config_reference)
    forced_lines = set(force_recs.keys())
    all_mismatches = []
    for d in wanted:
        driver_path = os.path.join(out_dirs[d], "input", "driver_main.txt")
        mism = d6_config_guard.check_against_reference(driver_path, config_reference, allow_lines=forced_lines)
        for lineno, ref_v, got_v in mism:
            all_mismatches.append((d, lineno, ref_v, got_v))
    if all_mismatches:
        if not args.allow_config_deviation:
            print("ERROR: D6 configuration guard -- MECHANISM record mismatch(es) against %s:" % config_reference)
            print(d6_config_guard.format_mismatches("*", [(l, r, g) for _, l, r, g in all_mismatches]))
            print("Refusing to run. Pass --allow-config-deviation to proceed (declared ablation arms only).")
            sys.exit(2)
        with open(os.path.join(out, "CONFIG_DEVIATION_LOG.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["deck", "record", "reference_value", "assembled_value"])
            w.writerows(all_mismatches)
        print("WARNING: D6 configuration guard -- %d MECHANISM-record deviation(s) logged to"
              " CONFIG_DEVIATION_LOG.csv (--allow-config-deviation set); proceeding." % len(all_mismatches))
    else:
        print("D6 configuration guard: all assembled decks match the v1.0 reference on every MECHANISM record.")

    print("%-20s %8s %8s %10s" % ("deck", "sec", "exit", "live_proc"))
    results = []
    any_fail = False
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(run_one, d, out_dirs[d], exe, exe_basename): d for d in wanted}
        for fut in as_completed(futs):
            d, dt, rc, n_live = fut.result()
            print("%-20s %8.1f %8d %10d" % (d, dt, rc, n_live))
            results.append((d, dt, rc, n_live))
            if rc != 0 or n_live > 0:
                any_fail = True

    md5_rows = []
    for d in wanted:
        ddir = out_dirs[d]
        md5_rows.append((d, "run.log", md5_file(os.path.join(ddir, "run.log"))))
        odir = os.path.join(ddir, "output")
        for root, _, files in os.walk(odir):
            for f in sorted(files):
                p = os.path.join(root, f)
                md5_rows.append((d, os.path.relpath(p, ddir), md5_file(p)))
    with open(os.path.join(out, "RUN_MD5.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["deck", "file", "md5"])
        w.writerows(md5_rows)

    if any_fail:
        print("FAIL: at least one deck returned non-zero exit or left a live process")
        sys.exit(1)
    print("All %d deck(s) completed, exit=0, 0 live processes." % len(wanted))


if __name__ == "__main__":
    main()
