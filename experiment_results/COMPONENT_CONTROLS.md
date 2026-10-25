# Four component controls: local full-schedule matching

This is the earlier Z-memory experiment, distinct from the common unitary-encoder 15 us benchmark. Original failure counts, primary DEMs and all recorded metadata are preserved.

| Arm | Failures / shots | LER (Wilson 95% CI) | Raw parity error | CZ | Padding | Full schedule |
|---|---:|---|---:|---:|---:|---:|
| 1 round + idle | 375/2000 | 18.75% [17.10%, 20.52%] | 33.45% | 19 | 8.24 us | 410.872 us |
| 1 round + CZ echo | 394/2000 | 19.70% [18.02%, 21.50%] | 35.50% | 67 | 7.904 us | 410.876 us |
| 1 round + ancilla measurement | 366/2000 | 18.30% [16.67%, 20.05%] | 34.45% | 19 | 6.168 us | 410.872 us |
| 3 rounds | 651/2000 | 32.55% [30.53%, 34.64%] | 41.70% | 67 | 0 us | 410.872 us |

| Contrast vs 1 round + idle | LER difference (pp) | Unadjusted 95% difference CI (pp) | Holm p (3 contrasts) |
|---|---:|---|---:|
| r1_cz_echo | +0.95 | [-1.49, +3.39] | 0.9403 |
| r1_ancilla_measure | -0.45 | [-2.86, +1.96] | 0.9403 |
| r3 | +13.80 | [+11.12, +16.46] | 4.406e-23 |

## Interpretation

- Three complete rounds have 13.80 percentage points higher LER than one round plus waiting; the gap has the same direction in both saved blocks.
- Extra paired-CZ echo and repeated ancilla measurement show no statistically detected difference here. This is not an equivalence test or evidence that gate/readout noise is absent.
- The CZ intervention is 48 extra gates in ideal CZ-CZ identity pairs on the 24 check links, not arbitrary unpaired entangling gates.
- The measurement intervention repeats the already-measured ancillas twice in Z, without re-entangling data or adding H. It records 16 auxiliary bits which are deliberately excluded from the primary decoder.
- All three one-round controls have the identical primary DEM; three rounds use their own full-history DEM.
- Full schedule durations, including initialization/readout, differ by at most 4 ns. Operation additions replace waiting; the contrast measures their net effect, not isolated error per gate.
- These interventions do not reproduce the full gate/measurement/basis-rotation sequence of repeated extraction, so the gap cannot be uniquely attributed to CZ or measurement alone.
- The original plan/results/job metadata also list an old unencoded r0 arm. Those records remain verbatim for provenance; its measurement data are not included and it is excluded from all submitted analyses.
- Missing isolated memory-boundary timings are blank in timing.csv; they are not inferred by subtracting an assumed reset time. Server timestamps and full waveforms remain unavailable.
