# Beam-Focusing for Physical Layer Security in Massive MIMO Systems

**Team:** Antoni, Jakub, Michał, Dawid, Filip

## About

Research code for studying Physical Layer Security (PLS) in Massive MIMO systems for 6G networks. Each scenario contrasts **sub-6 GHz Massive MIMO** (Nt = 32) with **mmWave Ultra-Massive MIMO** (Nt = 512) under passive and active eavesdropper attacks.

All simulations rely on MATLAB + Phased Array System Toolbox.

---

## Repository layout

```
.
├── run_all.m                   # one-click batch runner
├── simulations/                # one .m file per scenario
├── utils/                      # shared helpers (params, FSPL, ULA, ...)
└── results/                    # auto-generated white-background PNGs
```

### `utils/`
| File | Purpose |
|------|---------|
| `default_params.m`     | canonical project-wide parameters (carriers, antennas, SNR, RNG seed) |
| `setup_ula.m`          | half-wavelength ULA + phased.SteeringVector for a given Nt, fc |
| `compute_fspl.m`       | Free-Space Path Loss (Friis), linear and dB |
| `secrecy_rate.m`       | `max(0, R_bob - R_eve)` per Wyner |
| `jains_fairness.m`     | Jain's fairness index over per-user rates |
| `uniform_quantize.m`   | mid-rise uniform quantizer for complex baseband (low-resolution DAC model) |
| `jakes_correlation.m`  | Jakes temporal correlation J₀(2π f_d Δt) for Doppler-aged channels |
| `save_figure.m`        | export figure as 200 dpi PNG with white background to `/results` |

### `simulations/`
| Script | Topic | Key output |
|--------|-------|------------|
| `Ghz6_band_vs_mmWave_band.m`      | Baseline FSPL + MRT/ZF comparison @ 6 GHz vs 28 GHz | `fig_baseline_6GHz_vs_28GHz.png` |
| `sim_spatial_correlation.m`       | Vector vs Matrix ZF normalization under exponential correlation; Bob + Eve + Jain index | `fig_spatial_correlation.png` |
| `sim_colluding_eavesdroppers.m`   | Worst-case MRC attack by L cooperating Eves; 6 GHz vs 28 GHz | `fig_colluding_eavesdroppers.png` |
| `sim_artificial_noise.m`          | Null-space AN at 6 GHz; sweeps over L and over data-power fraction φ | `fig_artificial_noise.png` |
| `sim_pilot_contamination.m`       | Active beam-hijacking; 6 GHz vs 28 GHz, secrecy collapse + beam patterns | `fig_pilot_contamination.png` |
| `sim_fairness_normalization.m`    | Vector vs Matrix ZF: Sum-Rate vs Jain Fairness, sweeps over SNR and K | `fig_fairness_normalization.png` |
| `sim_phase_noise.m`               | Hardware phase-error robustness, 6 GHz vs 28 GHz, beam smearing | `fig_phase_noise.png` |
| `sim_channel_hardening.m`         | Channel-hardening proof: ‖h‖²/Nt → 1, Var ∝ 1/Nt, Secrecy Rate stability vs Nt | `fig_channel_hardening.png` |
| `sim_low_res_dac.m`               | Low-resolution DAC quantisation (b=1..5 bit), Bussgang factor, Bob/Eve/Secrecy vs bits | `fig_low_res_dac.png` |
| `sim_pilot_jamming.m`             | Denial-of-service jamming on the pilot phase; JPR sweep, perfect-CSI baseline, Nt heat-map | `fig_pilot_jamming.png` |
| `sim_csi_aging.m`                 | Doppler-aged CSI under Jakes model + L-tap Wiener predictor as defence | `fig_csi_aging.png` |
| `sim_location_error.m`            | Pointing-error sweep showing the narrow-beam paradox; 6 vs 28 GHz | `fig_location_error.png` |

---

## Running

```matlab
>> cd massive-mimo-pls-beamforming-main
>> run_all                                    % refresh every figure
% or run a single scenario:
>> addpath utils simulations
>> sim_spatial_correlation
```

All scripts use a fixed RNG seed (`default_params.rng_seed = 2026`) so results are reproducible.

---

## Conventions used across all scripts

* **SNR convention:** `SNR_dB` always denotes the *received* SNR after path loss; transmit power is derived as `P_tx = 10^(SNR_dB/10) * FSPL_lin`.
* **Noise:** `noise_var = 1` (linear), unless overridden in a scenario.
* **Secrecy Rate:** computed via the helper `secrecy_rate` (Wyner clamp at 0).
* **Plots:** every figure is saved to `/results` with white background via `save_figure`. Multi-panel layout (left = metric, right = physical insight) is the default.

---

## Roadmap

Sprint 1 ✅ done: Jain fairness for ZF normalization, hardware phase noise, channel hardening proof.
Sprint 2 ✅ done: Low-resolution DAC quantisation (Bussgang), pilot jamming (DoS attack).
Sprint 3 ✅ done: CSI aging (Jakes) + Wiener predictor; location-error amplification.
