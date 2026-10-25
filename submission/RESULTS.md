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

Timing values in tables/timing.csv are archived local pulse-compiler summaries, not server timestamps or full waveforms.
The native gate sequences match the direct fixed-layout no-routing construction and all CZ links exist in the saved architecture. This is a gate-level verification.
Near-zero/5 us/15 us results come from separate jobs. Calibration ID is shared; drift is not ruled out.
