# Config files

## Configs used for the current real-world analysis

The `_wml` configs (and `config_cdi_clean.yml`) turn on the nocovar RF pass
(`run_nocovar_for_rf: True`) and are the current versions.

| Dataset | Config |
|---|---|
| `autism_2_microbiomeHD` | `config_wml.yml` |
| `cdi_3_microbiomeHD_clean` | `config_cdi_clean.yml` |
| `ibd_3_CMD` | `config_IBD_wml.yml` |
| `crc_11_CMD` | `config_CRC_11_wml.yml` |

For `ibd_3_CMD` and `crc_11_CMD`, batch-corrected data and evaluation outputs exist in
more than one place (`trial_rw/`, `benchmark/benchmarked_results/`, `outputs/`). The
Oct 2026 percentile_norm rerun read both from `benchmark/benchmarked_results/`, while
`config_IBD_wml.yml` points to `trial_rw/benchmarked_results/`. Check
`post_integration_outputs` / `evaluation_outputs` before rerunning these two.

Simulation: `config_simulate_052626.yml`.

## Same filename, different contents

These files exist both at the repo root and in `additional_ymls/` with different contents:

- `config_cdi.yml`: both use `cdi_3_microbiomeHD`; the `additional_ymls/` copy has
  `condition_value: ASD`, which looks like a copy-paste error. The root copy (`CDI`) is the
  newer one.
- `config_CRC.yml`: root is `crc_8_CMD` (`adenoma`), `additional_ymls/` is `crc_3_CMD`
  (`CRC`). These are different datasets; neither is `crc_11_CMD`.
- `config_IBD.yml`: same dataset; the `additional_ymls/` copy is an older version with
  fewer methods and `./outputs` paths.

## Archived

`additional_ymls/archive/` holds `config_cdi_orig_rerun.yml` (one-off cdi rerun) and
`config_simulate_check.yml` (small test config).
