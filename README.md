# OTFS Channel Estimation

This is the code for [Beyond the Delay-Doppler Domain: A Time-Frequency Framework for Low-Overhead, Scalable OTFS Channel Estimation](https://arxiv.org/abs/2511.08504).

## Requirements

- MATLAB R2024b (verified) and above
- Communications Toolbox for QAM/PSK modulation, bit conversion, and AWGN
- Signal Processing Toolbox for peak-processing routines
- Parallel Computing Toolbox only when `use_parallel = true` and `trials > 1`

## Quick start

The package resolves its target pool relative to its own location, so the caller does not need to change the working directory.

```matlab
project_dir = '/path/to/Code_OTFS_ChannelEstimation_Modular';
addpath(project_dir);

params = config_params_std(struct( ...
    'SNR_db_list', 20, ...
    'trials', 1, ...
    'show_figures', false, ...
    'verbose_output', true));

results = run_sim_v2(params);
```

To run the default plotting script:

```matlab
run(fullfile(project_dir, 'run_single.m'));
```

The default full sweep is computationally expensive. Use one SNR point and one trial for an initial smoke test.

## Paper-specific code

`paper_2022_rg_bl_v2`, `paper_2023_ap_sip_bl`, and `paper_2024_online_bl_v2` contain the independently runnable simulations corresponding to their respective 2022, 2023, and 2024 studies.

The 2024 v2 configuration, initial RGBL estimator, online ORGBL estimator, and CSI metrics are included locally; no `paper_2024_online_bl` folder or sibling project is needed. Add the project root and `paper_2024_online_bl_v2` to the MATLAB path, then call `config_params_2024_v2(overrides)` and `run_sim_2024_v2(params)`. The 2024 online pilots are configured independently of the root pilot layout, including for small custom grids.

The 2024 defaults use `SNR_db_list = -10:5:30`, `trials = 2`, and `use_parallel = false`; caller overrides are preserved. Multi-trial runs skip detection and return `NaN` BER with `detection_skipped = true`. The root integer-delay channel and its target geometry are unchanged.

All three folders follow the root run controls below; their BER fields are `NaN` when detection is skipped. Each provides `run_single_*`, `run_experiment_*`, and `run_ber_batch_*`, and saves results in its own `results` folder.

The 2024 PHY supports unequal integer path delays without truncating shorter
path tails. To run its regression, add `paper_2024_online_bl_v2/tests` to the
MATLAB path and call `test_2024_multipath`.

## Configuration

`config_params_std` supports three target modes, in priority order:

1. `targets`: explicit deterministic target structs.
2. `n_targets`: random selection from the packaged pool.
3. `selected_targets`: fixed indices into the packaged pool.

Use `potential_targets` or `targets_file` to replace the packaged pool. Use `target_dopplers` or `doppler_range` to control target Doppler values.

Run controls, shared by `run_sim_v2` and the paper runners:

- `use_parallel`: run trials with `parfor` when `trials > 1`.
- `rng_seed_start`: trial `i` uses `rng(rng_seed_start + i - 1, 'threefry')` at every SNR, so serial and `parfor` runs draw identical realizations; `run_sim_v2` records the seeds in `results.rng_seeds`.
- `compute_ber`: empty runs detection only when `trials = 1`; `true` detects on every trial.

Multi-trial BER is expensive, so batched BER is a separate run: `run_ber_batch(overrides)` (optionally `run_ber_batch(overrides, @config_params_HF)`) sets `compute_ber = true` and saves to `results/ber_*.mat`. JSON experiments can set `"compute_ber": true` in `overrides`.

## Outputs

`run_sim_v2` returns:

- `NMSE_before`: channel NMSE after coarse estimation.
- `NMSE_after`: channel NMSE after successive refinement.
- `BER_after`: detected BER when detection runs (see `compute_ber`); `NaN` otherwise.
- `pct_err_coarse`: per-target coarse parameter errors.
- `pct_err_ssr`: per-target refined parameter errors.
- `n_detected_aoa`: detected coarse angle count per trial and SNR.
- `coarse_ests_all`: matched per-target coarse estimates.
- `SNR_db_list`: simulated SNR values.
- `params`: resolved configuration.
- `rng_seeds`: Monte Carlo seed used for each trial.

## Citation

```bibtex
@misc{wang2026delaydopplerdomaintimefrequencyframework,
  title        = {Beyond the Delay-Doppler Domain: A Time-Frequency Framework for Low-Overhead, Scalable OTFS Channel Estimation},
  author       = {Kailong Wang and Athina Petropulu},
  year         = {2026},
  eprint       = {2511.08504},
  archivePrefix = {arXiv},
  primaryClass = {eess.SP},
  url          = {https://arxiv.org/abs/2511.08504}
}
```

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).
