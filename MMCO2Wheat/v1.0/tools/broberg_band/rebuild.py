# Rebuild Broberg et al. 2019 Table 1 grain-yield response functions from their Supplementary Database
import pandas as pd, numpy as np
d = pd.read_excel("x/agronomy-497566-supplementary.xlsx", header=0)
d.columns = [str(c).strip() for c in d.columns]
print(d.shape); 
for c in []: #"Fumigation facility","CO2 treat.","O3 treat.","Temp. treat","Water treat.","Nitrogen treat.","Sowing time treat."]:
    print(c, d[c].value_counts(dropna=False).to_dict())
d = d[d["GY_mean"].notna()]
key = ["Reference","City","Year","Cultivar","Fumigation facility","Rooting env.","O3 treat.","Temp. treat","Water treat.","Sowing time treat.","Nitrogen treat."]
d[key] = d[key].astype(str)
pts = []
for k, g in d.groupby(key, dropna=False):
    if g["CO2 treat."].nunique() < 2: continue
    x, y = g["CO2 (ppm)"].astype(float).values, g["GY_mean"].astype(float).values
    b = np.polyfit(x, y, 1) if len(set(x)) > 1 else None
    if b is None: continue
    y350 = np.polyval(b, 350.0)
    for xi, yi in zip(x, y):
        pts.append((xi, 100 * (yi / y350 - 1)))
P = np.array(pts); print("points", len(P))
np.save("pts.npy", P)
for deg in (1, 2):
    print(deg, np.polyfit(P[:, 0], P[:, 1], deg)[::-1])
