# OTFS Channel Estimation

This is the code for [Beyond the Delay-Doppler Domain: A Time-Frequency Framework for Low-Overhead, Scalable OTFS Channel Estimation](https://arxiv.org/abs/2511.08504).

## Requirements

- MATLAB R2024b (verified) and above
- Communications Toolbox for QAM/PSK modulation, bit conversion, and AWGN
- Signal Processing Toolbox for peak-processing routines
- Parallel Computing Toolbox only when `trials > 1` and `params.n_workers > 0`

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

results = run_sim(params);
```

To run the default plotting script:

```matlab
run(fullfile(project_dir, 'run_single.m'));
```

The default full sweep is computationally expensive. Use one SNR point and one trial for an initial smoke test. `run_experiment('experiments/<name>.json')` runs a JSON-defined sweep.

## Paper-specific code

`paper_2022_rg_bl_v2`, `paper_2023_ap_sip_bl`, and `paper_2024_online_bl_v2` contain the independently runnable simulations corresponding to their respective 2022, 2023, and 2024 studies. Each provides `run_single_*`, `run_experiment_*`, and `run_ber_batch_*`, and saves results in its own `results` folder.

The 2024 v2 configuration, initial RGBL estimator, online ORGBL estimator, and CSI metrics are included locally; no `paper_2024_online_bl` folder or sibling project is needed. Add the project root and `paper_2024_online_bl_v2` to the MATLAB path, then call `config_params_2024_v2(overrides)` and `run_sim_2024_v2(params)`. The 2024 online pilots are configured independently of the root pilot layout, including for small custom grids.

The 2024 defaults use `SNR_db_list = -10:5:30`, `trials = 2`, and `use_parallel = false`; caller overrides are preserved. Multi-trial runs skip detection and return `NaN` BER with `detection_skipped = true`. The root integer-delay channel and its target geometry are unchanged.

The paper folders use their own run controls: `use_parallel` (run trials with `parfor` when `trials > 1`) and `compute_ber` (empty runs detection only when `trials = 1`; `true` detects on every trial).

The 2024 PHY supports unequal integer path delays without truncating shorter path tails.

## Configuration

`config_params_std` supports three target modes, in priority order:

1. `targets`: explicit deterministic target structs.
2. `n_targets`: random selection from the packaged pool.
3. `selected_targets`: fixed indices into the packaged pool.

Use `potential_targets` or `targets_file` to replace the packaged pool. Use `target_dopplers` or `doppler_range` to control target Doppler values. `include_los = true` adds a static direct path (known to the SSR, not detected by the coarse stage).

`angle_resolution`, `delay_resolution` and `doppler_resolution` round the true angles, delays and Dopplers and are the default SSR minimum steps (`ssr.min_step.a/r/v` override). Coarse estimates are on integer grids. `padLen` is the fixed guard length and must cover the largest delay.

Trial `i` uses `rng(rng_seed_start + i - 1, 'threefry')` at every SNR; trials run in parallel on `n_workers` workers when `trials > 1` (`0` = serial). BER is computed for `trials == 1`; `run_ber_batch('results/<file>.mat')` computes it for saved multi-trial results from the stored seeds and estimates.

## Outputs

`run_sim` returns:

- `err_coarse`, `err_ssr`: per-reflection `[angle, delay, Doppler]` errors, `(trials, SNR, target, 3)`; NaN = unmatched or failed.
- `nmse_coarse`, `nmse_ssr`: channel NMSE after coarse estimation and after refinement.
- `ests_coarse`, `ests_ssr`, `chanParams`: per-trial estimates and true channels.
- `n_detected_coarse`, `n_detected_ssr`, `n_matched_coarse`, `n_matched_ssr`: estimated and matched reflection counts.
- `BER`: detected BER for single-trial runs, else NaN.
- `rng_seeds`, `failed`, `SNR_db_list`, `params`.

Replot saved results with `plot_results(results)`.

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
