<!-- Copyright (c) 2018-2026 David Helman and The Hebrew University of Jerusalem -->
<!-- Licensed under the M&M-CO2Crop Research Licence v1.0 (see LICENSE): non-commercial research use with citation; redistribution of modified versions and commercial use require the copyright holders' consent. -->

# Known limitations of M&M-CO2Wheat v1.0 — optional diagnostics that are unavailable

M&M-CO2Wheat v1.0 ships one implemented formulation of each process, together with the
alternative formulations of Part 1 Table 6, which are present in the source and
documented as built-and-off. The source code of the other process options is not part of
this release, and a driver value that would have selected one is refused at startup
rather than accepted.

Some of those options were also the **writer of an optional diagnostic file**. For those,
the diagnostic file can never be produced by this release: the record that selected it is
refused, and the writer is not present. The affected file names, the records and the
released values are listed here, so that a user meets them in the documentation rather
than as a file that never appears.

None of this affects the 88 output files a v1.0 run writes, the physics, or any published
result. Every setting listed below is off in the released configuration.

## 1. Seven settings whose diagnostic writer is not present

In each row the record variable is the name the model prints when it refuses the value,
the released value is the only value the record accepts, and the last column gives the
unit-number declaration that survives in the source for the absent writer.

| # | driver record | record variable | unavailable setting | released value | diagnostic file that is unavailable | unit-number declaration in the source |
|---|---|---|---|---|---|---|
| 1 | 99 | `tiller_survival_mech` | mode 2 | 4 | `output/daily_tiller_afford.csv` | `main_MMCO2Crop.f90:388` (`unit_tca`) |
| 2 | 134 | `n_afford_priority_law` | 1 (and above) | 0 | `output/daily_n_afford_priority.csv` | `main_MMCO2Crop.f90:239` (`unit_nap`) |
| 3 | 135 | `n_payment_luxury_pass2` | 1 (and above) | 0 | `output/daily_n_payment_luxury_pass2.csv` | `main_MMCO2Crop.f90:241` (`unit_npl2`) |
| 4 | 136 | `n_afford_standing_basis` | 1 (and above) | 0 | `output/daily_n_afford_standing_basis.csv` | `main_MMCO2Crop.f90:243` (`unit_nsb`) |
| 5 | 137 | `floret_n_law` | 1 (and above) | 0 | `output/daily_floret_n.csv` | `main_MMCO2Crop.f90:245` (`unit_fnl`) |
| 6 | 139 | `n_uptake_diffusion` | mode 2 | 1 | `output/daily_n_uptake_bound.csv` | `main_MMCO2Crop.f90:251` (`unit_nub`) |
| 7 | 139 | `n_uptake_diffusion` | mode 3 | 1 | `output/daily_n_uptake_mode3.csv` | `main_MMCO2Crop.f90:253` (`unit_nuc`) |

Record 139 mode 1 — the released value — still writes `output/daily_n_uptake_terms.csv`,
which is one of the 88 files; only its mode-2 and mode-3 companions are unavailable.
Record 99 mode 4 — the released value — still writes
`output/daily_tiller_afford4.csv`, also one of the 88; only the mode-2 file is
unavailable.

Setting any of these seven values does not produce an empty file and does not fall
through to a different mechanism. The run is refused before the first simulated day, and
the model prints the record number, the record name, the offending value and the set of
implemented values, then exits non-zero.

## 2. Eight further conditional diagnostics that are also unavailable

The same situation, reached from records whose only accepted value is the one that leaves
the diagnostic switched off.

| driver record | record variable | released value | unavailable diagnostic file | `OPEN` statement in the source |
|---|---|---|---|---|
| 95 | `tleaf_class_diag` | 0 | `output/hourly_tleaf_classes.csv` | present, but unreachable at the released value |
| 102 | `stem_sheath_sep_mode` | 0 | `output/daily_stem_sep.csv` | absent |
| 106 | `wsc_law_member` | 0 | `output/daily_wsc_law.csv` | absent |
| 116 | `limb_census_mode` | 0 | `output/limb_census.csv` | absent |
| 132 | `leafn_hold_law` | 0 | `output/daily_leafn_hold.csv` | absent |
| 133 | `leafn_turnover_law` | 0 | `output/daily_leafn_turnover.csv` | absent |
| 152 | `census_files` | 1 | `output/ci_clamp_census.csv` | absent |
| 152 | `census_files` | 1 | `output/spac_census.csv` | absent |

Records 109-113 are the exception in this family: their Part 1 Table 6 alternative is
retained, so `output/internode_ranks.csv` and `output/internode_shortfall.csv` are
genuinely reachable with record 109 = 1, and both writers are present.

## 3. How to read this file

* A file named in this list is **not** a defect to report. Its absence is a documented
  consequence of releasing one configuration of each process plus the Part 1 Table 6
  alternatives.
* A diagnostic listed here cannot be recovered by editing the driver file. It requires
  the writer to be restored in the source.
* The accepted value of every restricted driver record is tabulated in section 2.2 of
  `../INPUT_GUIDE.md`.
