# apply.py <decks> <calibrated_phen.csv|calibrated_cult.csv>: write calibrated values to every deck of the group/cultivar
import sys, csv, os, importlib.util
spec = importlib.util.spec_from_file_location("fr", os.path.join(os.path.dirname(__file__), "face_run.py"))
os.environ.setdefault("FCAL_DECKS", sys.argv[1]); os.environ.setdefault("FCAL_EXE", "x")
fr = importlib.util.module_from_spec(spec); spec.loader.exec_module(fr)
from pathlib import Path
D = Path(sys.argv[1]); ALL = sorted(x for x in os.listdir(D) if (D / x).is_dir())
PHEN = {"WDF25_Ruta_aCO2_WW": "WDF25_Ruta_", "WDF25_Zahir_aCO2_WW": "WDF25_Zahir_", "WDF23_Ruta_aCO2": "WDF23_Ruta_",
        "WDF23_Zahir_aCO2": "WDF23_Zahir_",
        "BRAUNSCHWEIG_103": ["BRAUNSCHWEIG_%d" % i for i in (101, 102, 103, 104, 105, 106, 107, 110)],
        "BRAUNSCHWEIG_115": ["BRAUNSCHWEIG_%d" % i for i in (113, 114, 115, 116, 117, 118, 119, 122)],
        "MARICOPA_906": ["MARICOPA_%d" % i for i in (905, 906, 907, 908)], "MARICOPA_910": ["MARICOPA_910", "MARICOPA_912"],
        "MARICOPA_914": ["MARICOPA_%d" % i for i in (913, 914, 915, 916)]}
CULT = {"WDF25_Ruta_aCO2_WW": "_Ruta_", "WDF25_Zahir_aCO2_WW": "_Zahir_", "BRAUNSCHWEIG_103": "BRAUNSCHWEIG_", "MARICOPA_906": "MARICOPA_"}
for r in csv.DictReader(open(sys.argv[2])):
    c = r["deck"]
    if "Anthesis" in r:
        g = PHEN[c]; decks = [x for x in ALL if x.startswith(g)] if isinstance(g, str) else g
        r0 = fr.read_pheno(D / c / "input" / "phenology_cultivar.txt")
        new = fr.pheno_params(r0, float(r["Anthesis"]), float(r["GrainFill_TDU"]))
        for d in decks:
            assert abs(fr.read_pheno(D / d / "input" / "phenology_cultivar.txt")["Anthesis"] - r0["Anthesis"]) < 1e-6, d
            fr.set_rows(D / d / "input" / "phenology_cultivar.txt", new)
        print(c, "->", len(decks), "decks", {k: round(v, 1) for k, v in new.items()})
    else:
        decks = [x for x in ALL if CULT[c] in x or x.startswith(CULT[c])]
        for d in decks:
            if r.get("Blade_L_max_cm"):
                fr.set_rows(D / d / "input" / "phenology_cultivar.txt", {"Blade_L_max_cm": float(r["Blade_L_max_cm"])})
            if r.get("fruiting_efficiency") not in (None, "", "NA"):
                fr.set_fe(D / d / "input" / "cultivar_yield.csv", float(r["fruiting_efficiency"]))
        print(c, "->", len(decks), "decks", {k: r[k] for k in r if k not in ("deck", "runs")})
