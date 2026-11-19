#!/usr/bin/env python3
"""compare.py <produced_output_dir> <expected_output_dir> [--rtol R] [--atol A]

Compares every CSV file of an M&M-CO2Wheat output/ directory with the expected one.
Standard library only (Python 3.6 or newer).

Each numeric cell pair (a, b) is within tolerance when ANY of these holds:
  1. |a - b| <= atol + rtol * max(|a|, |b|)            (defaults rtol 1e-6, atol 1e-9)
  2. |a - b| <= 0.5 * 10**(-d), d = the number of decimals printed in the cell
     (for E format, the mantissa decimals minus the exponent); the smaller of the
     two cells' units is used
  3. the column is a named near-zero residual column (RESIDUAL_ATOL below) and
     |a - b| <= its absolute tolerance.
A cell outside tolerance in a named diagnostic column (DIAGNOSTIC below) is reported
as DIFFERENT-DIAGNOSTIC, with its file, column and row; it is not a FAIL.

File verdicts:
  IDENTICAL             the file is the same byte for byte
  PASS                  every cell is within tolerance
  DIFFERENT-DIAGNOSTIC  the only cells outside tolerance are in named diagnostic columns
  FAIL                  any other cell is outside tolerance, or rows, columns or text differ
  MISSING               the file was not produced
OVERALL is PASS only if no file is FAIL or MISSING.
"""
import csv, math, os, re, sys

# Near-zero residual columns: the difference of two large terms that is zero up to
# rounding. Absolute tolerance in the column's unit, with the reason.
RESIDUAL_ATOL = {
    # Soil heat closure G_surface_integral - storage_change. Its terms are typically
    # 1e5-1e6 J m-2 (daily median 4e5, hourly median 1.6e5), so double-precision
    # rounding leaves 1e-8 to 1e-6 J m-2. 1e-3 J m-2 is 1e-8 of the median daily term
    # and still flags any real closure error.
    ("daily_soil_heat.csv", "heat_residual_J_m2"): 1e-3,
    ("hourly_soil_heat.csv", "heat_residual_J_m2"): 1e-3,
}

# Diagnostic columns, checked in the source code to be written to output only and read
# by no state variable.
DIAGNOSTIC = {
    # c_actual = N paid / DM built that day (growth_module.f90). When the lamina increment
    # is a rounding remainder (about 1e-14 g m-2) instead of exactly 0, the ratio is the
    # target concentration instead of 0. It is read only by the output writer and by the
    # output-only counter lamina_below_target_cum.
    ("daily_n_construction.csv", "c_actual_L"): "N concentration of tissue built that day; output only",
    ("daily_n_construction.csv", "c_actual_S"): "as c_actual_L, stem; output only",
    ("daily_n_construction.csv", "c_actual_R"): "as c_actual_L, root; output only",
    # Difference of two values that are compared in their own columns (Tac/Tcan_demand_C
    # and Tair_C). It is stored only in the output buffer ebh_row and written once.
    ("hourly_energy_balance.csv", "Tac_minus_Tair_demand"): "Tac_demand_C - Tair_C; output only",
    ("hourly_energy_balance.csv", "Tcan_minus_Tair_demand"): "Tcan_demand_C - Tair_C; output only",
}

NUM = re.compile(r"^\s*[-+]?(\d*)\.?(\d*)(?:[eEdD]([-+]?\d+))?\s*$")


def num(s):
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except ValueError:
        return None


def unit(s):
    """Value of one unit in the last printed digit of the numeric string s."""
    m = NUM.match(s)
    if not m:
        return 0.0
    dec = len(m.group(2)) if "." in s else 0
    exp = int(m.group(3)) if m.group(3) else 0
    return 10.0 ** (exp - dec)


def compare(name, p, e, rtol, atol):
    """Return verdict, max_abs, max_rel, note, list of diagnostic cells."""
    if open(p, "rb").read() == open(e, "rb").read():
        return "IDENTICAL", 0.0, 0.0, "", []
    rp = list(csv.reader(open(p, newline="")))
    re_ = list(csv.reader(open(e, newline="")))
    if len(rp) != len(re_):
        return "FAIL", float("nan"), float("nan"), "row count %d vs %d" % (len(rp), len(re_)), []
    head = [h.strip() for h in re_[0]] if re_ else []
    mabs = mrel = 0.0
    bad = ""
    diag = []
    for i, (a_row, b_row) in enumerate(zip(rp, re_)):
        if len(a_row) != len(b_row):
            return "FAIL", mabs, mrel, "column count differs at row %d" % (i + 1), diag
        for j, (a, b) in enumerate(zip(a_row, b_row)):
            x, y = num(a), num(b)
            if x is None or y is None:
                if a.strip() != b.strip():
                    return "FAIL", mabs, mrel, "text differs at row %d col %d" % (i + 1, j + 1), diag
                continue
            d = abs(x - y)
            r = d / max(abs(x), abs(y), atol)
            if d > mabs: mabs = d
            if r > mrel: mrel = r
            col = head[j] if j < len(head) else str(j + 1)
            if d <= atol + rtol * max(abs(x), abs(y)):
                continue
            if d <= 0.5 * min(unit(a), unit(b)):
                continue
            if (name, col) in RESIDUAL_ATOL and d <= RESIDUAL_ATOL[(name, col)]:
                continue
            if (name, col) in DIAGNOSTIC:
                diag.append((col, i + 1, a.strip(), b.strip()))
                continue
            if not bad:
                bad = "first out-of-tolerance cell: row %d col %s (%s vs %s)" % (i + 1, col, a.strip(), b.strip())
    if bad:
        return "FAIL", mabs, mrel, bad, diag
    if diag:
        return "DIFFERENT-DIAGNOSTIC", mabs, mrel, "%d diagnostic cell(s), listed below" % len(diag), diag
    return "PASS", mabs, mrel, "", diag


def main():
    args = sys.argv[1:]
    rtol, atol = 1e-6, 1e-9
    if "--rtol" in args:
        k = args.index("--rtol"); rtol = float(args[k + 1]); del args[k:k + 2]
    if "--atol" in args:
        k = args.index("--atol"); atol = float(args[k + 1]); del args[k:k + 2]
    if len(args) != 2:
        print(__doc__); sys.exit(2)
    prod, exp = args
    files = sorted(f for f in os.listdir(exp) if f.endswith(".csv"))
    extra = sorted(set(f for f in os.listdir(prod) if f.endswith(".csv")) - set(files))
    counts = {}
    diags = []
    print("%-34s %-20s %12s %12s  %s" % ("file", "verdict", "max_abs", "max_rel", "note"))
    for f in files:
        p = os.path.join(prod, f)
        if not os.path.exists(p):
            v, a, r, note, dg = "MISSING", float("nan"), float("nan"), "not produced", []
        else:
            v, a, r, note, dg = compare(f, p, os.path.join(exp, f), rtol, atol)
        counts[v] = counts.get(v, 0) + 1
        diags += [(f,) + c for c in dg]
        print("%-34s %-20s %12.4g %12.4g  %s" % (f, v, a, r, note))
    for f in extra:
        print("%-34s %-20s %12s %12s  %s" % (f, "EXTRA", "", "", "not in the expected set (ignored)"))
    if diags:
        print()
        print("diagnostic cells outside tolerance (file, column, row: produced vs expected):")
        for f, col, row, a, b in diags:
            print("  %s, %s, row %d: %s vs %s" % (f, col, row, a, b))
    ok = counts.get("FAIL", 0) == 0 and counts.get("MISSING", 0) == 0
    print()
    print("files: %d expected; %s" % (len(files), ", ".join("%s %d" % kv for kv in sorted(counts.items()))))
    print("tolerance: rtol %g, atol %g; half a unit of the last printed digit; residual columns: %s" % (
        rtol, atol, ", ".join("%s %s %g" % (f, c, t) for (f, c), t in sorted(RESIDUAL_ATOL.items()))))
    print("OVERALL: %s" % ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
