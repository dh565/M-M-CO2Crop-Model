import csv,sys,os
R=sys.argv[1]
def col(f,name):
    r=list(csv.DictReader(open(f)))
    return [float(x[name]) for x in r if x[name].strip() not in ("","NaN")]
rows=[]
for co2 in (420,550):
  for tag in ("off","P40","P15","P8","P5","P3"):
    d=f"{R}/{tag}_{co2}/output"
    B=col(d+"/daily_fluxes_full.csv","Biomass_g_m2")[-1]
    Y=col(d+"/daily_yield_components.csv","grain_yield_g_m2")[-1]
    L=max(col(d+"/daily_d2_organs.csv","LAI_used"))
    U=pP=tp=cc=float("nan")
    if os.path.exists(d+"/daily_P.csv"):
        U=sum(col(d+"/daily_P.csv","uptake_gP_m2_d"))
        pP=min(col(d+"/daily_P.csv","p_paid_L")); tp=min(col(d+"/daily_P.csv","tpu_scale")); cc=min(col(d+"/daily_P.csv","c_cyt_mM"))
    rows.append((tag,co2,B,Y,L,U,pP,tp,cc))
print(f"{'case':5} {'CO2':>4} {'biomass':>8} {'grain':>7} {'LAImax':>6} {'Pupt kg/ha':>10} {'min pL':>6} {'minTPU':>6} {'min c_cyt':>9}")
for r in rows: print(f"{r[0]:5} {r[1]:4d} {r[2]:8.1f} {r[3]:7.1f} {r[4]:6.2f} {r[5]*10:10.2f} {r[6]:6.2f} {r[7]:6.2f} {r[8]:9.2f}")
print("\neCO2 gain (550 vs 420), % biomass / % grain:")
d={ (r[0],r[1]):r for r in rows}
for tag in ("off","P40","P15","P8","P5","P3"):
    a,e=d[(tag,420)],d[(tag,550)]
    print(f"{tag:5} {100*(e[2]/a[2]-1):6.1f} {100*(e[3]/a[3]-1):6.1f}")
