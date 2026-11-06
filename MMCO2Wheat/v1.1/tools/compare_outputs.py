#!/usr/bin/env python3
# Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem
# Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial
# research use with citation; redistribution of modified versions and commercial use
# require the copyright holders' consent.
"""
compare_outputs.py -- compare a produced M&M-CO2Wheat output/ directory against the
shipped expected output/ directory of the reference case.

Usage:
    python3 tools/compare_outputs.py <produced_dir> <expected_dir> [options]

Options:
    --csv <path>    also write the summary table as a CSV file
    --max-rows N    stop reading each file after N data rows (default: all)
    --quiet         print only the summary table and the final verdict

What it reports, per file:
    * byte-identical yes / no
    * if not byte-identical, for every numeric column: the maximum absolute
      difference, the maximum relative difference, and the data row (1-based,
      excluding the header) at which each maximum occurs
    * a one-line verdict

It then prints a summary table ordered by severity (worst maximum relative
difference first, maximum absolute difference as the tie-breaker) and a final
count of files in each verdict class.

Standard library only. No numpy, no pandas. Python 3.6+.

Verdict classes
    IDENTICAL   byte-for-byte equal
    NUMERIC     same shape, same non-numeric text, values differ only numerically
    STRUCTURE   row count, column count or header differs, or a text field differs
    MISSING     file present in expected/ but absent from produced/
    EXTRA       file present in produced/ but absent from expected/
    UNREADABLE  file could not be parsed as CSV on one of the two sides
"""

from __future__ import print_function

import csv
import math
import os
import sys

REL_FLOOR = 0.0  # relative difference uses |expected| as the denominator; see note below


def sorted_csv_names(d):
    out = []
    for name in sorted(os.listdir(d)):
        if name.startswith('.'):
            continue
        if os.path.isfile(os.path.join(d, name)):
            out.append(name)
    return out


def read_bytes(path):
    with open(path, 'rb') as fh:
        return fh.read()


def read_table(path, max_rows=None):
    """Return (header, rows). rows is a list of lists of strings."""
    with open(path, 'r', newline='') as fh:
        reader = csv.reader(fh)
        rows = []
        header = None
        for i, row in enumerate(reader):
            if i == 0:
                header = row
                continue
            rows.append(row)
            if max_rows is not None and len(rows) >= max_rows:
                break
    if header is None:
        header = []
    return header, rows


def as_float(tok):
    """Return float(tok), or None if tok is not a finite number."""
    t = tok.strip()
    if t == '':
        return None
    try:
        v = float(t)
    except ValueError:
        # Fortran sometimes writes 1.0D+03 or ***** on overflow
        t2 = t.replace('D', 'E').replace('d', 'e')
        try:
            v = float(t2)
        except ValueError:
            return None
    if math.isnan(v) or math.isinf(v):
        return None
    return v


class FileResult(object):
    def __init__(self, name):
        self.name = name
        self.verdict = None
        self.note = ''
        self.identical = False
        self.max_abs = 0.0
        self.max_abs_col = ''
        self.max_abs_row = 0
        self.max_rel = 0.0
        self.max_rel_col = ''
        self.max_rel_row = 0
        self.rel_undefined = 0   # cases where expected == 0 and produced != 0
        self.n_diff_cells = 0
        self.n_numeric_cols = 0
        self.columns = []        # per-column detail, worst first

    def sort_key(self):
        # worst first; UNREADABLE/STRUCTURE/MISSING/EXTRA rank above NUMERIC
        rank = {'MISSING': 5, 'EXTRA': 5, 'UNREADABLE': 5, 'STRUCTURE': 4,
                'NUMERIC': 3, 'IDENTICAL': 0}[self.verdict]
        return (-rank, -self.max_rel, -self.max_abs, self.name)


def compare_one(name, prod_path, exp_path, max_rows=None):
    r = FileResult(name)

    pb = read_bytes(prod_path)
    eb = read_bytes(exp_path)
    if pb == eb:
        r.verdict = 'IDENTICAL'
        r.identical = True
        r.note = 'byte-for-byte equal'
        return r

    try:
        ph, pr = read_table(prod_path, max_rows)
        eh, er = read_table(exp_path, max_rows)
    except Exception as exc:                      # noqa: BLE001
        r.verdict = 'UNREADABLE'
        r.note = 'CSV parse failed: %s' % (exc,)
        return r

    struct = []
    if ph != eh:
        if len(ph) != len(eh):
            struct.append('header column count %d vs expected %d' % (len(ph), len(eh)))
        else:
            struct.append('header text differs')
    if len(pr) != len(er):
        struct.append('data rows %d vs expected %d' % (len(pr), len(er)))

    ncol = min(len(ph), len(eh)) if (ph or eh) else 0
    nrow = min(len(pr), len(er))

    # per-column accumulators
    col_stats = {}
    for j in range(ncol):
        col_stats[j] = {'name': eh[j] if j < len(eh) else '(col %d)' % (j + 1),
                        'max_abs': 0.0, 'abs_row': 0,
                        'max_rel': 0.0, 'rel_row': 0,
                        'rel_undef': 0, 'n_diff': 0, 'numeric': False,
                        'text_diff': 0}

    for i in range(nrow):
        prow = pr[i]
        erow = er[i]
        if len(prow) != len(erow):
            struct.append('row %d field count %d vs expected %d'
                          % (i + 1, len(prow), len(erow)))
            continue
        ncell = min(len(prow), ncol) if ncol else len(prow)
        for j in range(ncell):
            pt = prow[j]
            et = erow[j]
            if pt == et:
                continue
            cs = col_stats.get(j)
            if cs is None:
                continue
            pv = as_float(pt)
            ev = as_float(et)
            if pv is None or ev is None:
                cs['text_diff'] += 1
                continue
            cs['numeric'] = True
            cs['n_diff'] += 1
            d = abs(pv - ev)
            if d > cs['max_abs']:
                cs['max_abs'] = d
                cs['abs_row'] = i + 1
            if ev != 0.0:
                rel = d / abs(ev)
                if rel > cs['max_rel']:
                    cs['max_rel'] = rel
                    cs['rel_row'] = i + 1
            elif d != 0.0:
                cs['rel_undef'] += 1

    cols = []
    for j in range(ncol):
        cs = col_stats[j]
        if cs['n_diff'] == 0 and cs['text_diff'] == 0:
            continue
        cols.append(cs)
        r.n_diff_cells += cs['n_diff']
        r.rel_undefined += cs['rel_undef']
        if cs['numeric']:
            r.n_numeric_cols += 1
        if cs['max_abs'] > r.max_abs:
            r.max_abs = cs['max_abs']
            r.max_abs_col = cs['name']
            r.max_abs_row = cs['abs_row']
        if cs['max_rel'] > r.max_rel:
            r.max_rel = cs['max_rel']
            r.max_rel_col = cs['name']
            r.max_rel_row = cs['rel_row']
    cols.sort(key=lambda c: (-c['max_rel'], -c['max_abs'], c['name']))
    r.columns = cols

    text_diffs = sum(c['text_diff'] for c in cols)
    if struct or text_diffs:
        r.verdict = 'STRUCTURE'
        parts = list(struct[:3])
        if text_diffs:
            parts.append('%d non-numeric field(s) differ' % text_diffs)
        r.note = '; '.join(parts)
    else:
        r.verdict = 'NUMERIC'
        r.note = ('%d cell(s) differ in %d column(s); max abs %.6g, max rel %.6g'
                  % (r.n_diff_cells, r.n_numeric_cols, r.max_abs, r.max_rel))
    return r


def fmt(x):
    if x == 0.0:
        return '0'
    return '%.4e' % x


def main(argv):
    args = [a for a in argv[1:]]
    opts = {'csv': None, 'max_rows': None, 'quiet': False}
    pos = []
    i = 0
    while i < len(args):
        a = args[i]
        if a == '--csv':
            i += 1
            opts['csv'] = args[i]
        elif a == '--max-rows':
            i += 1
            opts['max_rows'] = int(args[i])
        elif a == '--quiet':
            opts['quiet'] = True
        elif a in ('-h', '--help'):
            print(__doc__)
            return 0
        else:
            pos.append(a)
        i += 1

    if len(pos) != 2:
        sys.stderr.write('usage: compare_outputs.py <produced_dir> <expected_dir> '
                         '[--csv out.csv] [--max-rows N] [--quiet]\n')
        return 2

    prod_dir, exp_dir = pos[0], pos[1]
    for d in (prod_dir, exp_dir):
        if not os.path.isdir(d):
            sys.stderr.write('not a directory: %s\n' % d)
            return 2

    prod = sorted_csv_names(prod_dir)
    exp = sorted_csv_names(exp_dir)
    prod_set, exp_set = set(prod), set(exp)

    results = []
    for name in exp:
        if name not in prod_set:
            r = FileResult(name)
            r.verdict = 'MISSING'
            r.note = 'not produced'
            results.append(r)
            continue
        results.append(compare_one(name,
                                   os.path.join(prod_dir, name),
                                   os.path.join(exp_dir, name),
                                   opts['max_rows']))
    for name in prod:
        if name not in exp_set:
            r = FileResult(name)
            r.verdict = 'EXTRA'
            r.note = 'produced but not expected'
            results.append(r)

    print('M&M-CO2Wheat v1.0 -- output comparison')
    print('produced : %s' % os.path.abspath(prod_dir))
    print('expected : %s' % os.path.abspath(exp_dir))
    print('files in expected/ : %d        files in produced/ : %d' % (len(exp), len(prod)))
    print('')

    if not opts['quiet']:
        print('--- per-file detail (differing files only) -------------------------------')
        any_detail = False
        for r in sorted(results, key=lambda x: x.sort_key()):
            if r.verdict == 'IDENTICAL':
                continue
            any_detail = True
            print('')
            print('%s  [%s]' % (r.name, r.verdict))
            print('  byte-identical : no')
            print('  %s' % r.note)
            if r.rel_undefined:
                print('  %d cell(s) where the expected value is exactly 0 and the produced '
                      'value is not: relative difference undefined, see max abs'
                      % r.rel_undefined)
            if r.columns:
                print('  %-34s %12s %8s %12s %8s' %
                      ('column', 'max_abs', 'at_row', 'max_rel', 'at_row'))
                for c in r.columns:
                    print('  %-34s %12s %8d %12s %8d' %
                          (c['name'][:34], fmt(c['max_abs']), c['abs_row'],
                           fmt(c['max_rel']), c['rel_row']))
        if not any_detail:
            print('  (none -- every expected file is byte-identical)')
        print('')

    print('--- summary, worst first -------------------------------------------------')
    print('%-40s %-11s %12s %12s %10s' %
          ('file', 'verdict', 'max_abs', 'max_rel', 'cells_diff'))
    for r in sorted(results, key=lambda x: x.sort_key()):
        print('%-40s %-11s %12s %12s %10d' %
              (r.name[:40], r.verdict, fmt(r.max_abs), fmt(r.max_rel), r.n_diff_cells))

    counts = {}
    for r in results:
        counts[r.verdict] = counts.get(r.verdict, 0) + 1
    print('')
    print('--- verdict counts -------------------------------------------------------')
    for k in ('IDENTICAL', 'NUMERIC', 'STRUCTURE', 'MISSING', 'EXTRA', 'UNREADABLE'):
        if k in counts:
            print('  %-11s %4d' % (k, counts[k]))
    n_ident = counts.get('IDENTICAL', 0)
    print('')
    print('%d of %d expected file(s) byte-identical.' % (n_ident, len(exp)))
    if n_ident == len(exp) and not counts.get('EXTRA', 0):
        print('VERDICT: byte-identical reproduction of the reference case.')
        rc = 0
    elif counts.get('MISSING') or counts.get('STRUCTURE') or counts.get('UNREADABLE'):
        print('VERDICT: NOT reproduced -- structural differences present. Read the '
              'per-file detail above before interpreting any number.')
        rc = 2
    else:
        worst = max((r.max_rel for r in results), default=0.0)
        worst_abs = max((r.max_abs for r in results), default=0.0)
        print('VERDICT: numerically different, same structure. Worst maximum relative '
              'difference %s, worst maximum absolute difference %s.'
              % (fmt(worst), fmt(worst_abs)))
        print('         This is the expected situation on a different compiler, '
              'compiler version, optimization level or platform. No measured '
              'cross-platform tolerance exists for this model; see docs/START_HERE.md.')
        rc = 1

    if opts['csv']:
        with open(opts['csv'], 'w', newline='') as fh:
            w = csv.writer(fh)
            w.writerow(['file', 'verdict', 'byte_identical', 'max_abs',
                        'max_abs_column', 'max_abs_row', 'max_rel',
                        'max_rel_column', 'max_rel_row', 'cells_differing',
                        'rel_undefined_cells', 'note'])
            for r in sorted(results, key=lambda x: x.sort_key()):
                w.writerow([r.name, r.verdict, 'yes' if r.identical else 'no',
                            repr(r.max_abs), r.max_abs_col, r.max_abs_row,
                            repr(r.max_rel), r.max_rel_col, r.max_rel_row,
                            r.n_diff_cells, r.rel_undefined, r.note])
        print('')
        print('summary CSV written: %s' % os.path.abspath(opts['csv']))

    return rc


if __name__ == '__main__':
    sys.exit(main(sys.argv))
