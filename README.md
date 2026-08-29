# OTFS Channel Estimation

Self-contained MATLAB implementation of the primary OTFS channel-estimation, sensing-assisted refinement, and data-detection workflow from `OTFS_CE_Modula_v2`. Paper-specific reproduction implementations are intentionally excluded from the public package.

> **Maintainer space:** Replace this note with a one-paragraph project description before release.

## Contents

| Path | Purpose |
|---|---|
| `run_single.m` | Default simulation, result display, and plots |
| `run_sim_v2.m` | SNR-by-trial Monte Carlo workflow used by `run_single.m` |
| `config_params_std.m` | Default 4 GHz system, pilot, target, and SSR configuration |
| `targets/` | Packaged target pools |

All runtime source and default target data for the primary workflow are included. Directories matching `paper_*/`, `OTFS_sample_code/`, `experiments/Temp/`, generated results, local tooling state, archives, and MATLAB scratch files are excluded by `.gitignore`.

## Requirements

- MATLAB R2024b (verified)
- Communications Toolbox for QAM/PSK modulation, bit conversion, and AWGN
- Signal Processing Toolbox for peak-processing routines
- Parallel Computing Toolbox only when `trials > 1`


## Quick start

The package resolves bundled target files relative to its own location, so the caller does not need to change the working directory.

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
BER is evaluated when `trials == 1`. Multi-trial runs intentionally skip the expensive detector and return zero in `BER_after` as the sentinel for “not computed”; do not average those sentinel values as measured BER.



## Reproducibility

`run_sim_v2` seeds trial `i` with `rng_seed_start + i - 1`. The returned `results.rng_seeds` records the sequence. Example:

```matlab
params = config_params_std(struct( ...
    'SNR_db_list', 20, ...
    'trials', 1, ...
    'rng_seed_start', 101));
results = run_sim_v2(params);
```

## Release metadata

- **License:** BSD 3-Clause; see [`LICENSE`](LICENSE).
- **Citation:** [Add preferred paper citations and BibTeX.]
- **Authors/contact:** [Add public maintainer details.]
- **Repository:** <https://github.com/kwang0913/Code_OTFS_ChannelEstimation_Modular>
- **Validation:** MATLAB R2024b on Linux; compact serial and two-trial parallel smoke workflows.

## Citation

[Add the preferred citation text and BibTeX here.]

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).
