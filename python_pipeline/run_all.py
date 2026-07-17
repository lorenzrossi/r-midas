#!/usr/bin/env python3
# =============================================================================
# PYTHON PIPELINE - ORCHESTRATOR (equivalent of run_all.R)
#
# DEFAULT (current project focus): AR benchmark + updated standard R-MIDAS
# (Almon over the AR lag set J = {1,2,3,7}), then the summary table and the
# combination / subperiod / GR evaluation restricted to the r_midas family.
#
# The other stages (r_midas_ext, ru_midas, surveys_h) remain registered and
# can be selected explicitly, e.g.:
#   STAGES=ar,r_midas,ru_midas,summary python3 run_all.py
#
# Env knobs (forwarded to every stage): COUNTRY_NAME, SPEC_LIMIT, N_CORES,
# WINDOW_DAYS, RMIDAS_DISCOUNT, ORIGIN_LIMIT, FAMILIES.
# =============================================================================

import importlib
import os
import time

STAGE_MODULES = {
    "ar": "ar_benchmark",
    "r_midas": "r_midas",
    "r_midas_ext": "r_midas_extended",
    "ru_midas": "ru_midas",
    "surveys_h": "r_midas_surveys_h",
    "summary": "summary_table",
    "combination": "combination_subperiod_gr",
}


def _load_stage(stage):
    """Import only the requested stage, so optional files may be absent."""
    module = importlib.import_module(STAGE_MODULES[stage])
    return module.main


# Current focus: AR + updated R-MIDAS only (other stages via STAGES env var).
DEFAULT_ORDER = ["ar",
                 "r_midas",
                 "summary",
                 "combination"]


def main():
    stages = os.environ.get("STAGES", "")
    order = ([s.strip() for s in stages.split(",") if s.strip()]
             if stages else DEFAULT_ORDER)
    for st in order:
        if st not in STAGE_MODULES:
            print(f"Unknown stage '{st}' (known: {list(STAGE_MODULES)})")
            continue
        print(f"\n{'='*70}\nSTAGE: {st}\n{'='*70}")
        # Restrict the combination stage to the r_midas family unless the
        # user set FAMILIES explicitly.
        old_families = os.environ.get("FAMILIES")
        set_families_here = st == "combination" and old_families is None
        if set_families_here:
            os.environ["FAMILIES"] = "r_midas"
        t0 = time.time()
        try:
            _load_stage(st)()
        finally:
            if set_families_here:
                os.environ.pop("FAMILIES", None)
        print(f"stage {st} done in {time.time() - t0:.1f}s")
    print("\nAll requested stages done.")


if __name__ == "__main__":
    main()
