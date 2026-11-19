# Table 4 rows and CO2 / water / nitrogen effects

These scripts produce the Part 2 Table 4 rows that `make_part2_tables.py` does not write: the treatment-level series rows and the extra rows (stage dates and the remaining harvest rows). They also produce the treatment-ratio effects quoted in the text.

## Steps

1. Run all 38 decks with the v1.0 executable:
   `bash run_decks.sh <mmco2crop> ../../decks <runs dir> [workers]`
2. Set the environment variables and run the scripts:
   ```
   export REF_OUTPUTS=<runs dir> OBS_DIR=<observations dir> OUT_DIR=<output dir>
   python table4_series.py ..  table4_series.json
   python table4_extra_rows.py ..  table4_extra_rows.json
   python effects_CO2_water_N.py ..  effects.md v1.0=<runs dir>
   ```

Here `..` is the `paper_scripts` folder; the scripts import `config.py` and `lib/` from it.
