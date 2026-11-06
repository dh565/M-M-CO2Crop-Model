# P module tests (v1.1, private branch v1.1-P-module)

- `run_co2xP.sh <run_dir>`: runs the reference case (WDF25_Ruta_aCO2_WW) with P off and with the P module on at soil-test P 40, 15, 8, 5, 3 mg kg-1, each at 420 and 550 umol mol-1 CO2.
- `summarize_co2xP.py <run_dir>`: table of biomass, grain, LAI max, P uptake, lowest paid fraction, lowest TPU scale, and the eCO2 gain per P level.

Run outputs go OUTSIDE the repository (scratch), never into ~/Documents or ~/Desktop.

## Result, 6 Nov 2026 (first code version, provisional parameters)

| soil-test P | eCO2 gain, biomass | eCO2 gain, grain | biomass at 420, g m-2 |
|---|---|---|---|
| P off | 9.6 % | 12.3 % | 1461 |
| 40 | 9.6 % | 12.3 % | 1461 (identical to off) |
| 15 | 9.5 % | 12.3 % | 1460 |
| 8 | 7.8 % | 10.1 % | 1430 |
| 5 | -1.2 % | -6.2 % | 302 |
| 3 | 0.2 % | 0.0 % | 154 |

The loss of the CO2 gain under P shortage emerges; no CO2 x P term is coded.
Mass balances: soil < 1e-13, plant < 1e-15 g P m-2.
Against the field data of Rodriguez et al. (2000, Bray-1 5.5 mg kg-1: -70 % biomass, -73 % LAI at 61 d after emergence), the model at 5 mg kg-1 gives -44 % biomass and -31 % LAI at 61 d and 0 % at 44 d. The effect comes too late, because the tiller and phyllochron routes (phase B) are not coded yet.

## Revision, 6 Nov 2026 (later): measured tissue P composition, leaf appearance and tiller buds

- The ribosome cost from the MAXIMUM protein-synthesis rate (Raven 2013) gave about 10x too little RNA against the measured nucleic-acid P (Veneklaas 2012 p.309, 0.3-2 mg P g-1 DW; Sulpice 2014 ribosome counts). It is replaced by the measured composition of new tissue (nucleic acids, lipid, esters, cytosolic Pi). RNA decays to 25 % as tissue matures. The organ cost scales with the organ's protoplasm N (one N:P).
- Leaf appearance runs at the lamina P-paid fraction; tiller buds grow out in proportion to the P paid in their window.
- P40 = P off (no regression). Grain loading no longer drains the leaf cytosol.
- The supply side now decides everything. eps_test (the soil-test recovery of sorbed P) is unmeasured:
  - eps = 1.0: P15 collapses. Wrong: Bray-1 15 is near-sufficient for wheat.
  - eps = 0.3: P15 and P8 are sufficient; P5 gives -28 % biomass at 61 DAE (field: -70 %, Rodriguez 2000). P40 then takes up 93 kg P/ha, an unrealistic luxury uptake: transporters are repressed only by cytosolic Pi, which stays at its set point while the vacuole fills.
- Next: (1) vacuole-full feedback on transporter synthesis (whole-cell Pi signalling); (2) a measured eps for Bray-1 (isotopic E-value studies); (3) topsoil P stratification (the soil test is 0-20 cm); (4) mycorrhizal uptake (Li et al. 2006: >50 % of wheat P).

## Revision: tissue-Pi feedback on transporter synthesis

- Pht1 synthesis is now repressed by the tissue Pi (cytosol + vacuole). It was repressed by the clamped cytosolic Pi before.
- The replete reference is 5.6 mM, from the Veneklaas 2012 pool ratio.
- eps = 0.3, P40: uptake 93 -> 57 kg P ha-1 (plant P about 3.9 mg g-1 at maturity, high end of measured); P15 and P8 are unchanged (sufficient); P5: -2 % at 61 DAE, -41 % at 80 DAE.
- eps = 1.0: almost unchanged (P40 22.6 kg ha-1).
- The CO2 gain is unchanged where P is sufficient (9.6 %) and lost where P limits.
- The timing at moderately low P still depends on the unmeasured eps.
