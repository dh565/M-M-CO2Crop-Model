# Table 4 rows and CO2 / water / nitrogen effects

These scripts produce the Part 2 Table 4 rows that `make_part2_tables.py` does not write: the treatment-level series rows, the leaf gas-exchange and thermal-image rows, and the extra rows (stage dates and the remaining harvest rows). They also produce numbers quoted in the text: the treatment-ratio effects, the CO2 x water and CO2 x nitrogen interactions of the observations, the share of grain-filling days on which the protein-coupled bound sets grain Fe loading, and the soil mineral N at maturity at Braunschweig without and with the crop-driven immobilization.

## Steps

1. Run all 38 decks with the v1.0 executable:
   `bash run_decks.sh <mmco2crop> ../../decks <runs dir> [workers]`
2. Set the environment variables and run the scripts:
   ```
   export REF_OUTPUTS=<runs dir> OBS_DIR=<observations dir> OUT_DIR=<output dir>
   python table4_series.py ..  table4_series.json
   python table4_extra_rows.py ..  table4_extra_rows.json
   python effects_CO2_water_N.py ..  effects.md v1.0=<runs dir>
   python table4_agg.py ..  table4_agg.json
   python anova_summary.py ..  anova_summary.txt
   python fe_bound_days.py ..  fe_bound_days.txt
   python soil_nmin_immobilization.py ..  <mmco2crop> <work dir> soil_nmin_immobilization.txt [workers]
   ```

| script | output | where it is quoted |
|---|---|---|
| `table4_series.py` | treatment-level series rows (energy balance, soil water and temperature) | Table 4 |
| `table4_agg.py` | leaf An, gs, Tr rows of every site and the WheatDryFACE 2024/25 thermal-image leaf temperature row; pooled per-date n and r; per-date CO2 effects | Table 4; Sect. 6.6.1 and Fig. 3g-l |
| `table4_extra_rows.py` | stage dates and the remaining harvest rows | Table 4 |
| `effects_CO2_water_N.py` | CO2, water and nitrogen effects as treatment ratios, observed and simulated | Part 2 text (CO2, water and nitrogen effects) |
| `anova_summary.py` | F and P of the CO2 x water and CO2 x nitrogen interaction terms of the observations | Sect. 6.7 |
| `fe_bound_days.py` | days on which the protein-coupled bound sets grain Fe loading, over the days with grain growth, by site | Sects. 6.5 and 7.6 |
| `soil_nmin_immobilization.py` | soil mineral N 0-40 cm at maturity, eight Braunschweig decks, observed and simulated without (driver record 197 = 0, run by the script) and with the immobilization (the runs in `REF_OUTPUTS`, decks as shipped); RMSD of each | Sect. 4 (the fitted immobilization constant) |

`soil_nmin_immobilization.py` runs the eight decks itself: it copies each from `../../decks`, sets record 197 to 0 and runs the executable in `<work dir>/<deck>`.

Here `..` is the `paper_scripts` folder; the scripts import `config.py` and `lib/` from it.
