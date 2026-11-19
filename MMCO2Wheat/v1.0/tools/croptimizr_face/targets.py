# targets.py <ref_run> <out.json>: observed anthesis/maturity dates (control decks), grain yield and peak LAI (FE/leaf decks)
import sys, os, io, json, contextlib
sys.path.insert(0, os.getcwd()); sys.path.insert(0, 'lib')
import config
from lib import paperfigs as pf
with contextlib.redirect_stdout(io.StringIO()): pf.init(sys.argv[1])
P2 = pf.P2
CTRL = ["WDF25_Ruta_aCO2_WW","WDF25_Zahir_aCO2_WW","WDF23_Ruta_aCO2","WDF23_Zahir_aCO2","BRAUNSCHWEIG_103","BRAUNSCHWEIG_115","MARICOPA_906","MARICOPA_910","MARICOPA_914"]
out = {"phen": {}, "cult": {}}
for d in CTRL:
    st = {s: od for s, md, od, diff, scored, src in P2.PHEN[d]}
    out["phen"][d] = {"Anthesis": st.get("Anthesis"), "Maturity": st.get("Maturity")}
for d in ["WDF25_Ruta_aCO2_WW","WDF25_Zahir_aCO2_WW","BRAUNSCHWEIG_103","MARICOPA_906"]:
    g = lambda q: (P2.OBS_REG.get(d, q) or {}).get("value")
    out["cult"][d] = {"grain_yield": g("grain_yield"), "LAI_peak": g("LAI_peak"), "model_grain_yield": P2.MODELS[d].v.get("grain_yield")}
print(json.dumps(out, indent=1)); json.dump(out, open(sys.argv[2], "w"), indent=1)
