# Simulation scenarios

This document describes every script under `simulations/`: channel model, precoding, what is being studied, and the figure each script writes to `results/`.

See also: [CURATED_RESULTS.md](CURATED_RESULTS.md).

---

## Shared setup


| Item                | Convention                                                                                                   |
| ------------------- | ------------------------------------------------------------------------------------------------------------ |
| **Secrecy rate**    | Per-user R_s = \max(0, R_b - R_e); sum where noted (`secrecy_rate.m`, Wyner clamp)                           |
| **SNR**             | `SNR_rx_dB` is received power at Bob; abstract scripts use normalized Rayleigh (no FSPL in the channel draw) |
| **Noise**           | `noise_var = 1` unless a script overrides it                                                                 |
| **Reproducibility** | `default_params.rng_seed = 2026`                                                                             |
| **Startup**         | `pls_startup()` clears workspace, checks dependencies, closes figures                                        |
| **Requirements**    | MATLAB R2016b+; **Phased Array System Toolbox** (`assert_requirements`, `list_project_requirements`)         |
| **Figures**         | Dark-theme PNG via `save_figure` / `apply_plot_style` (Bob = green, Eve = salmon)                            |


### Channel model families


| Family                          | Implementation                                             | Used in                                 |
| ------------------------------- | ---------------------------------------------------------- | --------------------------------------- |
| **3GPP TR 38.901–inspired CDL** | `channel_3gpp_ula` + ULA steering (`setup_ula`)            | 7 dual-band or 6 GHz scenarios          |
| **i.i.d. Rayleigh**             | `(randn + 1j*randn)/sqrt(2)`                               | Fairness, hardening, DAC, pilot jamming |
| **Correlated Rayleigh**         | Exponential spatial correlation matrix R_{ij}=\rho^{|i-j|} | Spatial correlation                     |
| **Jakes temporal fading**       | `jakes_correlation` + Cholesky synthesis                   | CSI aging                               |


**CDL profiles** (not a full 3GPP system simulator):

- **6 GHz (`sub6`)** — CDL-A-like clusters, Rician K \approx 7 dB, richer NLOS.
- **28 GHz (`mmwave`)** — CDL-D-like clusters, K \approx 18 dB, stronger LOS.

---

## Scenario reference

### 1. `Ghz6_band_vs_mmWave_band` → `fig_baseline_6GHz_vs_28GHz.png`


|                       |                                                                                                                |
| --------------------- | -------------------------------------------------------------------------------------------------------------- |
| **Channel**           | 3GPP CDL-A @ 6 GHz (Nt = 32), CDL-D @ 28 GHz (Nt = 512); angles fixed (Bob / Eve)                              |
| **Precoding**         | MRT toward Bob; ZF with null space toward Eve                                                                  |
| **Sweep**             | Received SNR at Bob                                                                                            |
| **Research question** | How does the 6 GHz vs mmWave array size and cluster profile change secrecy under the same received-SNR budget? |
| **What to look for**  | mmWave narrow spatial structure vs sub-6 spatial diversity; MRT vs ZF gap; beam-pattern snapshots              |


---

### 2. `sim_moving_bob` → `fig_moving_bob.png`


|                       |                                                                                                                 |
| --------------------- | --------------------------------------------------------------------------------------------------------------- |
| **Channel**           | 3GPP CDL, both bands; Bob angle swept, Eve at fixed off-boresight angle                                         |
| **Precoding**         | **Fixed** beam steered to boresight (no re-tracking); MRT-style fixed beam vs ZF null toward Eve                |
| **Sweep**             | Bob bearing \theta_b                                                                                            |
| **Research question** | When the BS does not re-steer, how fast does secrecy collapse if Bob leaves the main lobe—especially at 28 GHz? |
| **What to look for**  | Sharp secrecy cliff for mmWave vs wider 6 GHz footprint; beam plot with Bob markers                             |


---

### 3. `sim_pilot_contamination` → `fig_pilot_contamination.png`


|                       |                                                                                                                |
| --------------------- | -------------------------------------------------------------------------------------------------------------- |
| **Channel**           | 3GPP CDL, 6 GHz vs 28 GHz                                                                                      |
| **Attack**            | Eve contaminates channel estimate: \hat{h} = h_b + \sqrt{\beta} h_e + n_{\mathrm{est}}; beam w \propto \hat{h} |
| **Sweep**             | Contamination strength \beta \in [0,1]                                                                         |
| **Research question** | Active pilot contamination (“beam hijacking”): when does Eve pull the beam toward herself and kill secrecy?    |
| **What to look for**  | Secrecy vs \beta; clean vs hijacked beam patterns; stronger sensitivity at 28 GHz                              |


---

### 4. `sim_colluding_eavesdroppers` → `fig_colluding_eavesdroppers.png`


|                       |                                                                                     |
| --------------------- | ----------------------------------------------------------------------------------- |
| **Channel**           | 3GPP CDL, 6 GHz vs 28 GHz; K Bobs, L Eves at random angles                          |
| **Precoding**         | ZF across Bobs; Eve side uses colluding MRC (sum of                                 |
| **Sweep**             | Number of cooperating Eves L                                                        |
| **Research question** | Worst-case colluding eavesdroppers: does ultra-massive mmWave still win as L grows? |
| **What to look for**  | Secrecy sum-rate and Jain fairness on secrecy rates vs L                            |


---

### 5. `sim_artificial_noise` → `fig_artificial_noise.png`


|                       |                                                                                                        |
| --------------------- | ------------------------------------------------------------------------------------------------------ |
| **Channel**           | 3GPP CDL-A @ 6 GHz only; K users, L Eves                                                               |
| **Precoding**         | ZF data in H null-space; AN in null(H^T) with power split \phi (data) / 1-\phi (AN)                    |
| **Sweep**             | (A) L at fixed \phi; (B) \phi at fixed L                                                               |
| **Research question** | Null-space artificial noise (Goh & Hong): trade-off between data power and Eve jamming under collusion |
| **What to look for**  | Optimal \phi; AN gain vs no-AN baseline                                                                |


---

### 6. `sim_location_error` → `fig_location_error.png`


|                       |                                                                                                                             |
| --------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| **Channel**           | 3GPP CDL, 6 GHz vs 28 GHz                                                                                                   |
| **Precoding**         | Beam from **estimated** angle \hat{\theta}*b = \theta_b + \mathcal{N}(0,\sigma*{\mathrm{loc}}^2)                            |
| **Sweep**             | Localization error std \sigma_{\mathrm{loc}}                                                                                |
| **Research question** | Narrow-beam paradox: mmWave needs accurate pointing; how does mis-pointing affect secrecy and “beam on target” probability? |
| **What to look for**  | Secrecy drop vs \sigma_{\mathrm{loc}}; 28 GHz loses target lock faster than 6 GHz                                           |


---

### 7. `sim_phase_noise` → `fig_phase_noise.png`


|                       |                                                                                        |
| --------------------- | -------------------------------------------------------------------------------------- |
| **Channel**           | 3GPP CDL, 6 GHz vs 28 GHz                                                              |
| **Impairment**        | W_{\mathrm{err}} = W \odot e^{j\sigma_\phi \mathcal{N}(0,I)} on precoder columns       |
| **Precoding**         | Matrix vs vector ZF normalization                                                      |
| **Sweep**             | Phase-error std (degrees)                                                              |
| **Research question** | Per-RF-chain phase jitter smears narrow beams; is mmWave more vulnerable at large N_t? |
| **What to look for**  | Secrecy vs \sigma_\phi; beam-pattern degradation snapshots                             |


---

### 8. `sim_spatial_correlation` → `fig_spatial_correlation.png`


|                       |                                                                                                                               |
| --------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| **Channel**           | Correlated Rayleigh: H = R^{1/2} H_{\mathrm{iid}}, R_{ij}=\rho^{|i-j|}; 12 dB user gain spread; attenuated Eve                |
| **Precoding**         | Matrix (Frobenius) vs vector (per-user) ZF normalization                                                                      |
| **Sweep**             | Spatial correlation \rho                                                                                                      |
| **Research question** | Antenna correlation breaks symmetry between the two ZF power constraints—effect on secrecy and **Jain fairness on Bob rates** |
| **What to look for**  | Secrecy sum-rate and fairness vs \rho; vector norm usually fairer                                                             |


---

### 9. `sim_fairness_normalization` → `fig_fairness_normalization.png`


|                       |                                                                                                                                                  |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Channel**           | i.i.d. Rayleigh + 12 dB gain spread; attenuated Eve (`eve_attn_dB`)                                                                              |
| **Precoding**         | Matrix vs vector ZF (see [CURATED_RESULTS.md](CURATED_RESULTS.md#zf-normalization--fairness-sim_fairness_normalization-sim_spatial_correlation)) |
| **Sweep**             | SNR and number of users K                                                                                                                        |
| **Research question** | Sum-rate vs fairness trade-off between global Frobenius budget and per-stream power cap                                                          |
| **What to look for**  | Jain index on **Bob** rates (not secrecy); vector normalization typically higher Jain                                                            |


---

### 10. `sim_channel_hardening` → `fig_channel_hardening.png`


|                       |                                                                                                                                |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| **Channel**           | i.i.d. Rayleigh h \sim \mathcal{CN}(0,I_{N_t}) for Bob and Eve                                                                 |
| **Precoding**         | MRT / matched filtering style secrecy evaluation                                                                               |
| **Sweep**             | N_t                                                                                                                            |
| **Research question** | Channel hardening: |h|^2/N_t \to 1 reduces Bob’s rate variance; does Eve harden equally, and what happens to secrecy variance? |
| **What to look for**  | Histograms at small vs large N_t; \mathrm{Var}(|h|^2/N_t) and std of secrecy vs N_t                                            |


---

### 11. `sim_low_res_dac` → `fig_low_res_dac.png`


|                       |                                                                                                        |
| --------------------- | ------------------------------------------------------------------------------------------------------ |
| **Channel**           | i.i.d. Rayleigh multi-user + Eve                                                                       |
| **Impairment**        | Uniform DAC quantization (`uniform_quantize`); Bussgang factor for effective SNR                       |
| **Sweep**             | DAC bits (1–5, ideal), SNR                                                                             |
| **Research question** | Low-resolution DACs (6G hardware): how many bits are needed to preserve secrecy? (Xu et al., TWC 2018) |
| **What to look for**  | Bob / Eve / secrecy vs bits; empirical vs theoretical Bussgang factor                                  |


---

### 12. `sim_pilot_jamming` → `fig_pilot_jamming.png`


|                       |                                                                                    |
| --------------------- | ---------------------------------------------------------------------------------- |
| **Channel**           | i.i.d. Rayleigh; noisy LS estimate \hat{H} = H + E + J                             |
| **Attack**            | Gaussian jamming on **training** slot (not pilot replay)                           |
| **Sweep**             | Jamming-to-pilot ratio (JPR); optional N_t heat-map                                |
| **Research question** | Training-phase DoS: how does pilot jamming degrade ZF secrecy vs perfect CSI?      |
| **What to look for**  | Secrecy collapse vs JPR; contrast with `sim_pilot_contamination` (hijack vs noise) |


---

### 13. `sim_csi_aging` → `fig_csi_aging.png`


|                       |                                                                                                                                  |
| --------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| **Channel**           | Jakes-correlated fading over time slots; 6 GHz, N_t = 32                                                                         |
| **CSI**               | (1) stale h[0]; (2) L-tap Wiener predictor; (3) perfect h[\tau] oracle                                                           |
| **Sweep**             | UE velocity; predictor order L                                                                                                   |
| **Research question** | Doppler-aged CSI hurts BS beamforming but not Eve’s local estimate—can Wiener prediction recover secrecy? (Zhu et al., TWC 2018) |
| **What to look for**  | Secrecy vs speed; gap between no prediction, Wiener, and perfect CSI                                                             |


---

## Topology maps

`generate_topologies.m` draws link-geometry sketches (not simulation results) into `topology/` — one PNG per scenario. Requires the same toolbox check as scenarios (`physconst` via `default_params`).

---

## Run order

`run_all.m` executes all 13 scripts above and writes `results/run_all.log`. Order matches the table in [README.md](README.md#simulations).