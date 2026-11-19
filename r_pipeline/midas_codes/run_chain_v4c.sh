#!/bin/bash
# =============================================================================
# TODO_v4c run chain (November 2026).  Sequential; every step logs into its
# run folder.  Launch from midas_codes/:
#   caffeinate -i -s bash run_chain_v4c.sh 2>&1 | tee chain_v4c_$(date +%Y%m%d_%H%M).log
# Steps:
#   1. density stage with DENSITY_RECAL=1 (R = 365) into results_recal/, and the
#      sensitivity R = 730 into results_recal730/  (Germany: estimation_v4_crps,
#      Italy: estimation_v4_it_log)                                       item 1
#   2. combination / summary / tests / cssed reruns of the existing runs, so the
#      combinations average the 252 macro specs (COMB_MACRO_ONLY=1)      item 4a
#   3. compare_vst_runs.R for Italy (tests csv)                            item 4c
#   4. check: R combinations vs comb252.py                                 item 4a
#   5. Germany full plain-asinh run (252 specs) + compare_vst_runs.R       item 3
#   6. paper tables (tables_v3/scripts), SPA / MCS in the foreground       item 5
# =============================================================================
MC="/Users/lorenz/Documents/Claude/Projects/Forecasting electricity prices/r-midas/r_pipeline/midas_codes"
PS="/Users/lorenz/Documents/Claude/Projects/Forecasting electricity prices/tables_v3/scripts"
export N_CORES=${N_CORES:-9}
T0=$(date +%s)
step()   { echo; echo "[$(date '+%F %T')] === $* ==="; }
runall() { local d="$1" lg="$2"; shift 2
           ( cd "$MC/$d" && env "$@" Rscript run_all.R > "$lg" 2>&1 ); local rc=$?
           echo "  [$(date '+%H:%M')] $d ($lg): exit $rc"; grep -E "stage '.*' (done|FAILED)|pipeline done" "$MC/$d/$lg" | sed 's/^/    /'; }
py()     { local lg="log_$(echo "$*" | tr ' ' '_').txt"; echo "-- python3 $*"; ( cd "$PS" && python3 "$@" > "$lg" 2>&1 ); echo "   exit $? ($lg)"; tail -2 "$PS/$lg" | sed 's/^/    /'; }

step "1. density stage: rescaled predictive sd (R=365) and sensitivity (R=730)"
runall estimation_v4_crps   recal_de.log    COUNTRY_NAME=Germany VST_MODE=normalized DENSITY_RECAL=1 RESULT_SUBDIR=results_recal    STAGES=density
runall estimation_v4_it_log recal_it.log    COUNTRY_NAME=Italy   VST_MODE=log        DENSITY_RECAL=1 RESULT_SUBDIR=results_recal    STAGES=density
runall estimation_v4_crps   recal730_de.log COUNTRY_NAME=Germany VST_MODE=normalized DENSITY_RECAL=1 RESULT_SUBDIR=results_recal730 RECAL_WINDOW=730 STAGES=density
runall estimation_v4_it_log recal730_it.log COUNTRY_NAME=Italy   VST_MODE=log        DENSITY_RECAL=1 RESULT_SUBDIR=results_recal730 RECAL_WINDOW=730 STAGES=density

step "2. combinations over the 252 macro specs: rerun combination, summary, tests, cssed"
FAM=arx,r_midas,r_midas_l7
runall estimation_v4          rerun252.log COUNTRY_NAME=ALL   VST_MODE=normalized FAMILIES=$FAM STAGES=combination,summary,tests,cssed
runall estimation_v4_crps     rerun252.log COUNTRY_NAME=ALL   VST_MODE=normalized STAGES=tests
runall estimation_v4_it_log   rerun252.log COUNTRY_NAME=Italy VST_MODE=log        FAMILIES=$FAM STAGES=combination,summary,tests,cssed
runall estimation_v4_it_level rerun252.log COUNTRY_NAME=Italy VST_MODE=level      FAMILIES=$FAM STAGES=combination,summary,tests,cssed
runall estimation_v4_it_asinh rerun252.log COUNTRY_NAME=Italy VST_MODE=asinh      FAMILIES=$FAM STAGES=combination,summary,tests,cssed

step "3. compare the Italian transformations (metrics + tests csv)"
( cd "$MC" && COMPARE_CC=it VST_BASE=normalized Rscript compare_vst_runs.R > compare_vst_it_v4c.log 2>&1 ); echo "  exit $?"
grep -E "^wrote|^Error" "$MC/compare_vst_it_v4c.log" | sed 's/^/    /'; ls -la "$MC/results_compare/" | sed 's/^/    /'

step "4. check 4a: R combinations (252 specs) vs comb252.py"
py comb252.py
py check_comb252.py

step "5. Germany: full plain-asinh run on the 252 macro specs"
mkdir -p "$MC/estimation_v4_de_asinh"
( cd "$MC/estimation_v4_de_asinh" && for f in ../*.R; do ln -sf "$f" .; done && ln -sfn ../data_v4 data )
runall estimation_v4_de_asinh run_de_asinh.log COUNTRY_NAME=Germany VST_MODE=asinh SPEC_REQUIRE_LF=1
( cd "$MC" && COMPARE_CC=de VST_BASE=normalized VST_RUNS="normalized=estimation_v4,asinh=estimation_v4_de_asinh" Rscript compare_vst_runs.R > compare_vst_de.log 2>&1 ); echo "  exit $?"
grep -E "^wrote|^Error" "$MC/compare_vst_de.log" | sed 's/^/    /'

step "6. paper tables (tables_v3/scripts)"
py comb252.py
py tab_main.py
py tab_win.py
py tab_marg.py
py tab_crps.py
py fig_gr.py
py gr_t.py
py spa.py de 1 7 14 21 28
py spa.py it 1 7 14 21 28
py mcs.py de
py mcs.py it
py mcs_summary.py
py tab_spa.py
py tab_app.py
py tab_vst.py

echo; echo "[$(date '+%F %T')] === CHAIN DONE in $(( ($(date +%s) - T0) / 60 )) min ==="
