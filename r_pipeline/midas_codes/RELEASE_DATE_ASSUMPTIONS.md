# Release-date assumptions for the low-frequency variables

The pipeline maps each monthly observation to the date on which it would
normally have become public, then carries the latest released value forward
day by day (`make_daily_lf` in `common_utils.R` / `common_utils.py`). Exact
historical release dates take precedence when the monthly dataset contains a
column named `<variable>_release_date` (e.g. `Manufacturing_release_date`,
`pmi_mfg_de_release_date`).

## Default conventional schedules

| variable | assumed release | env override (default) |
|---|---|---|
| Manufacturing PMI (`pmi_mfg_de`, `pmi_mfg_it`) | 1st working day of the following month | `LF_PMI_MFG_WORKDAY` (1) |
| Services PMI (`pmi_serv_de`, `pmi_serv_it`) | 3rd working day of the following month | `LF_PMI_SERV_WORKDAY` (3) |
| ifo Business Climate (`ifo_bci_de`) | 25th of the reference month | `LF_IFO_DAY` (25) |
| ISTAT business confidence (`istat_bci_it`) | 27th of the reference month | `LF_ISTAT_BCI_DAY` (27) |
| IPI components (`Consumer_Goods`, `Manufacturing`, `Energy`) | 45 calendar days after reference-month end | `LF_IPI_DELAY_DAYS` (45) |

Dates falling on weekends are rolled forward to the next Monday. These rules
avoid anticipatory use of a value but remain a pseudo-real-time
approximation when vintage-specific release dates are unavailable.

## Compatibility mode

Set `LF_RELEASE_MODE=calendar_lag` to restore the former fixed-month-shift
mapping (each day of month m uses the value of month m − pub_lag, with
pub_lag = 1 for surveys and 2 for IPI).
