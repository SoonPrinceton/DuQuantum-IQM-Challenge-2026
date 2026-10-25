# Reproduced measured results

All counts below are recomputed from the archived shots using their original DEM and MWPM settings. No tuning or new experiments.

| Dataset | Arm | Failures / shots | LER | Wilson 95% CI |
|---|---|---:|---:|---|
| memory_15us | Z_r0 | 502/2000 | 25.10% | [23.25%, 27.05%] |
| memory_15us | Z_r1 | 646/2000 | 32.30% | [30.29%, 34.38%] |
| memory_15us | Z_r2 | 730/2000 | 36.50% | [34.42%, 38.63%] |
| memory_15us | Z_r3 | 748/2000 | 37.40% | [35.31%, 39.54%] |
| memory_15us | Z_single | 13/2000 | 0.65% | [0.38%, 1.11%] |
| memory_15us | X_r0 | 970/2000 | 48.50% | [46.31%, 50.69%] |
| memory_15us | X_r1 | 909/2000 | 45.45% | [43.28%, 47.64%] |
| memory_15us | X_r2 | 998/2000 | 49.90% | [47.71%, 52.09%] |
| memory_15us | X_r3 | 972/2000 | 48.60% | [46.41%, 50.79%] |
| memory_15us | X_single | 591/2000 | 29.55% | [27.59%, 31.59%] |
| memory_near_zero | Z_r0_near0 | 58/2000 | 2.90% | [2.25%, 3.73%] |
| memory_near_zero | X_r0_near0 | 40/2000 | 2.00% | [1.47%, 2.71%] |
| memory_5us | Z_r0_5us | 234/2000 | 11.70% | [10.36%, 13.18%] |
| memory_5us | X_r0_5us | 1045/2000 | 52.25% | [50.06%, 54.43%] |
| component_controls | r1_idle | 375/2000 | 18.75% | [17.10%, 20.52%] |
| component_controls | r1_cz_echo | 394/2000 | 19.70% | [18.02%, 21.50%] |
| component_controls | r1_ancilla_measure | 366/2000 | 18.30% | [16.67%, 20.05%] |
| component_controls | r3 | 651/2000 | 32.55% | [30.53%, 34.64%] |

Timing values in tables/timing.csv are archived local pulse-compiler summaries, not server timestamps or full waveforms.
The primary 14 native gate sequences match the direct no-routing recipe; the four earlier controls reproduce their saved fixed-patch intervention recipe. All CZ links exist in the saved architecture. This is gate-level verification.
The component controls match full local schedules (~410.872 us), not the 15 us memory window. See COMPONENT_CONTROLS.md for their separate timing definition and contrasts.
Near-zero/5 us/15 us results come from separate jobs. Calibration ID is shared; drift is not ruled out.
