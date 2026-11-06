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
