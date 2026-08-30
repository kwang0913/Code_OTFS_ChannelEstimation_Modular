# OTFS Channel Estimation

This is the code for [Beyond the Delay-Doppler Domain: A Time-Frequency Framework for Low-Overhead, Scalable OTFS Channel Estimation](https://arxiv.org/abs/2511.08504).

## Requirements

- MATLAB R2024b (verified) and above
- Communications Toolbox for QAM/PSK modulation, bit conversion, and AWGN
- Signal Processing Toolbox for peak-processing routines
- Parallel Computing Toolbox only when `trials > 1`

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

## Configuration

`config_params_std` supports three target modes, in priority order:

1. `targets`: explicit deterministic target structs.
2. `n_targets`: random selection from the packaged pool.
3. `selected_targets`: fixed indices into the packaged pool.

Use `potential_targets` or `targets_file` to replace the packaged pool. Use `target_dopplers` or `doppler_range` to control target Doppler values.

`run_sim_v2` seeds trial `i` with `rng_seed_start + i - 1` and records the used seeds in `results.rng_seeds`.

## Outputs

`run_sim_v2` returns:

- `NMSE_before`: channel NMSE after coarse estimation.
- `NMSE_after`: channel NMSE after successive refinement.
- `BER_after`: detected BER for single-trial runs; zero is the “not computed” sentinel for multi-trial runs.
- `pct_err_coarse`: per-target coarse parameter errors.
- `pct_err_ssr`: per-target refined parameter errors.
- `n_detected_aoa`: detected coarse angle count per trial and SNR.
- `coarse_ests_all`: matched per-target coarse estimates.
- `SNR_db_list`: simulated SNR values.
- `params`: resolved configuration.
- `rng_seeds`: Monte Carlo seed used for each trial.

Do not average multi-trial `BER_after` zero sentinels as measured BER.

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
