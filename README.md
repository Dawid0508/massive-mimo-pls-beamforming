# Beam-Focusing for Physical Layer Security in Massive MIMO Systems

**Team:** Antoni, Jakub, Michał, Dawid, Filip

## About

Research code for studying Physical Layer Security (PLS) in Massive MIMO systems for 6G networks. Scenarios compare **sub-6 GHz Massive MIMO** (Nt = 32) with **mmWave Ultra-Massive MIMO** (Nt = 512) under passive and active eavesdropper attacks, plus abstract models for hardware, fairness, and CSI effects.

All simulations use **MATLAB R2016b or later** with the add-on below.

### Required software

| Component | Used for |
|-----------|----------|
| **MATLAB** (R2016b+) | Base linear algebra, `besselj`, plotting; `strings` in `run_all` |
| **Phased Array System Toolbox** | `phased.ULA`, `phased.SteeringVector`, `physconst`, `step(sv,...)` in `setup_ula` |
| **5G Toolbox** | `nrCDLChannel` — CDL-A channel realizations in all scenario scripts |

Install add-ons: **Home → Add-Ons → Get Add-Ons** → search *Phased Array System Toolbox* and *5G Toolbox* → Install → restart MATLAB.

Checks run from `run_all`, `pls_startup`, `default_params`, `setup_ula`, and `generate_topologies`. Dependency report:

```matlab
addpath utils
list_project_requirements
```

---

## Repository layout

```
.
├── run_all.m                   # batch all scenarios; log → results/run_all.log
├── generate_topologies.m       # link-geometry sketches → topology/
├── SCENARIOS.md                # per-scenario models, research goals, outputs
├── CURATED_RESULTS.md          # recommended figures for reports
├── simulations/                # one .m file per scenario
├── utils/                      # shared helpers (params, FSPL, ULA, 3GPP CDL, ...)
├── results/                    # simulation figures (dark-theme PNG, 200 dpi)
└── topology/                   # per-scenario top-down geometry maps
```

### Documentation

| File | Contents |
|------|----------|
| [SCENARIOS.md](SCENARIOS.md) | All 13 scenarios: channel model, precoding, sweeps, research questions |
| [CURATED_RESULTS.md](CURATED_RESULTS.md) | Which PNGs to highlight in a report |

### `utils/` (selected)

| File | Purpose |
|------|---------|
| `default_params.m` | Project-wide parameters (carriers, antennas, SNR, RNG seed) |
| `setup_ula.m` | Half-wavelength ULA + `phased.SteeringVector` |
| `channel_3gpp_ula.m` | TR 38.901-inspired CDL-A / CDL-D on ULA |
| `assert_requirements.m` | Fail fast if MATLAB or Phased Array Toolbox is missing |
| `list_project_requirements.m` | Print dependency audit |
| `secrecy_rate.m` | `max(0, R_bob - R_eve)` (Wyner) |
| `apply_plot_style.m` | Dark theme, line colors, legend + reference-line labels |
| `pls_colors.m` | Bob/Eve/band palette for plots |
| `save_figure.m` | Export dark-theme PNG via `apply_plot_style` |
| `mark_bob.m` / `mark_eve.m` | Bob/Eve xline markers (styled on save) |
| `pls_startup.m` | Per-scenario init + requirement check |

### `simulations/`

| Script | Output figure |
|--------|---------------|
| `Ghz6_band_vs_mmWave_band.m` | `fig_baseline_6GHz_vs_28GHz_nrCDL.png` |
| `sim_moving_bob.m` | `fig_moving_bob_doppler.png` |
| `sim_pilot_contamination.m` | `fig_pilot_contamination.png` |
| `sim_colluding_eavesdroppers.m` | `fig_colluding_eavesdroppers.png` |
| `sim_artificial_noise.m` | `fig_artificial_noise_tradeoffs.png` |
| `sim_location_error.m` | `fig_location_error.png` |
| `sim_phase_noise.m` | `fig_phase_noise.png` |
| `sim_spatial_correlation.m` | `fig_spatial_mrt_zf.png` |
| `sim_fairness_normalization.m` | `fig_fairness_normalization_cdl.png` |
| `sim_channel_hardening.m` | `fig_channel_hardening_cdl.png` |
| `sim_low_res_dac.m` | `fig_low_res_dac.png` |
| `sim_pilot_jamming.m` | `fig_pilot_jamming.png` |
| `sim_csi_aging.m` | `fig_csi_aging.png` |

Details for each row: [SCENARIOS.md](SCENARIOS.md).

---

## Running

```matlab
run_all                 % all scenarios → results/
generate_topologies     % geometry maps → topology/

% single scenario:
addpath utils simulations
sim_spatial_correlation
```

All scripts use `default_params.rng_seed = 2026` for reproducibility.

---

## Conventions

- **SNR:** `SNR_rx_dB` is received SNR at Bob (`rx_snr_power`). `print_scenario_snr` logs TX/RX per actor when `dist_m` / `fc_Hz` are set. Abstract Rayleigh scripts use normalized channels (no FSPL in H).
- **3GPP channels:** `channel_3gpp_ula` — CDL-A-like @ 6 GHz, CDL-D-like @ 28 GHz (cluster + ULA steering; not a full 3GPP stack).
- **Noise:** `noise_var = 1` (linear) unless overridden.
- **Secrecy rate:** `max(0, R_b - R_e)` per user; sum where noted.
- **Plots:** Dark background (`#1a1a2e`), light axes; Bob green, Eve salmon. `results/run_all.log` is git-ignored; PNGs are committed.

---

## Development history

- Sprint 1: Jain fairness / ZF normalization, phase noise, channel hardening.
- Sprint 2: Low-resolution DAC (Bussgang), pilot jamming.
- Sprint 3: CSI aging (Jakes + Wiener predictor), location-error / narrow-beam effects.
