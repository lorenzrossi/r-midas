# estimation_v4_de_asinh — Germany only, plain asinh transformation (TODO_v4c, item 3)

Same v4 data (`data -> ../data_v4`) and settings as `estimation_v4` (theta2 free, 7 starts, window
3650, all families, all stages), except `VST_MODE=asinh` (plain `asinh(P)`, the v3 transform) instead
of the window-normalised asinh, and `SPEC_REQUIRE_LF=1`: only the 252 specifications with at least one
monthly indicator are estimated (item 0; the daily-only alternatives are the AR models with gas / Brent).
Combinations therefore average the same 252 specs as the rerun combinations of `estimation_v4`.

Purpose: make the appendix comparison of transformations symmetric with Italy (`tables_v3/tab_app_vst.tex`,
Germany panel).  Compare with `estimation_v4` on the EUR scale only:
`results_compare/compare_vst_de.csv` and `compare_vst_tests_de.csv` (`compare_vst_runs.R`,
`VST_RUNS="normalized=estimation_v4,asinh=estimation_v4_de_asinh"`).

Launched by `../run_chain_v4c.sh`; log: `run_de_asinh.log`.
