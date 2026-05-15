# Recommended figures

Full scenario descriptions: [SCENARIOS.md](SCENARIOS.md).

Run all simulations: `>> run_all`  
Run geometry sketches only: `>> generate_topologies`

---

## Primary results (`results/`)


| Figure                            | Scenario            | Main takeaway                                              |
| --------------------------------- | ------------------- | ---------------------------------------------------------- |
| `fig_baseline_6GHz_vs_28GHz.png`  | Baseline            | 6 vs 28 GHz, CDL-A/D, MRT vs ZF, SNR sweep                 |
| `fig_moving_bob.png`              | Moving Bob          | Fixed beam; mmWave secrecy cliff when Bob leaves boresight |
| `fig_pilot_contamination.png`     | Pilot contamination | Active hijack: secrecy vs Eve pilot strength               |
| `fig_location_error.png`          | Location error      | Pointing error hurts narrow mmWave beams most              |
| `fig_colluding_eavesdroppers.png` | Colluding Eves      | Secrecy vs L cooperating eavesdroppers                     |
| `fig_artificial_noise.png`        | Artificial noise    | Null-space AN: data/noise split \phi at 6 GHz              |


---

## Additional results


| Figure                           | Topic                                                 |
| -------------------------------- | ----------------------------------------------------- |
| `fig_phase_noise.png`            | Phase noise; matrix vs vector ZF                      |
| `fig_channel_hardening.png`      | |h|^2/N_t concentration vs N_t                        |
| `fig_spatial_correlation.png`    | Correlation \rho; ZF normalization + Jain (Bob rates) |
| `fig_fairness_normalization.png` | SNR and K sweeps; sum-rate vs Jain (Bob rates)        |
| `fig_low_res_dac.png`            | DAC resolution vs secrecy                             |
| `fig_pilot_jamming.png`          | Training-slot jamming (DoS)                           |
| `fig_csi_aging.png`              | Doppler + Wiener CSI predictor                        |


---

## Topology maps (`topology/`)

Top-down geometry per scenario (BS, users, Eve, optional beam wedge):  
`topology_baseline.png` … `topology_moving_bob.png` (13 files).

---

## ZF normalization notes

Used in `sim_fairness_normalization` and `sim_spatial_correlation`:


| Mode       | Constraint            | Typical effect                     |
| ---------- | --------------------- | ---------------------------------- |
| **Matrix** | |W|*F^2 = P*{rx}      | Can boost weak ZF columns          |
| **Vector** | |W(:,k)|^2 = P_{rx}/K | More equal Bob rates → higher Jain |


- **Jain index** is computed on **Bob** rates R_b, not secrecy rates.
- **12 dB** user gain spread separates the two normalizations under Rayleigh fading.
- Eve is attenuated by `default_params.eve_attn_dB` (15 dB) in abstract runs.

Helpers: `zf_precoder_normalize`, `compute_zf_secrecy_metrics`, `attenuate_eve_channel`.

---

## SNR and channels

- **Logging:** `print_scenario_snr` — SNR_tx and SNR_rx per actor when geometry is set.
- **Abstract Rayleigh:** no FSPL in the channel; received SNR equals configured linear power.
- **3GPP-style:** `channel_3gpp_ula.m` — CDL-A @ 6 GHz, CDL-D @ 28 GHz.

