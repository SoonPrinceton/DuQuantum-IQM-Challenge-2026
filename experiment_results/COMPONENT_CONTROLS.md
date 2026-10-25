# Four component controls: first-syndrome-record window matching

Z-memory, distinct from the common unitary-encoder 15 us benchmark. Measured counts, primary DEMs and acquisition metadata are preserved.
First complete syndrome record to final readout barrier, target 5.348 us; first three arms rerun, r3 is historical; NOT a common independent encoding boundary

| Arm | Failures / shots | LER (Wilson 95% CI) | Raw parity error | CZ | Padding | Full schedule |
|---|---:|---|---:|---:|---:|---:|
| 1 round + idle | 279/2000 | 13.95% [12.50%, 15.54%] | 32.30% | 19 | 5.348 us | 407.98 us |
| 1 round + CZ echo | 285/2000 | 14.25% [12.79%, 15.85%] | 30.85% | 67 | 5.008 us | 407.976 us |
| 1 round + ancilla measurement | 250/2000 | 12.50% [11.12%, 14.02%] | 29.50% | 19 | 3.276 us | 407.98 us |
| 3 rounds | 651/2000 | 32.55% [30.53%, 34.64%] | 41.70% | 67 | 0 us | 410.872 us |

| Contrast vs 1 round + idle | LER difference (pp) | Unadjusted 95% difference CI (pp) | Holm p (3 contrasts) |
|---|---:|---|---:|
| r1_cz_echo | +0.30 | [-1.86, +2.46] | 0.8203 |
| r1_ancilla_measure | -1.45 | [-3.55, +0.65] | 0.3824 |
| r3 | +18.60 | [+16.03, +21.14] | 2.661e-44 |

## Interpretation

- Three-round minus one-round-idle LER difference: +18.60 percentage points. Refer to the intervals and adjusted tests above.
- A non-significant contrast is not an equivalence test or evidence that gate/readout noise is absent.
- The CZ intervention is 48 extra gates in ideal CZ-CZ identity pairs on the 24 check links, not arbitrary unpaired entangling gates.
- The measurement intervention repeats the already-measured ancillas twice in Z, without re-entangling data or adding H. It records 16 auxiliary bits which are deliberately excluded from the primary decoder.
- All three one-round controls have the identical primary DEM; three rounds use their own full-history DEM.
- Matching uses the first complete syndrome record, not a common independently encoded-state barrier. The reused three-round measurements come from earlier jobs; drift is a confound.
- These interventions do not reproduce the full gate/measurement/basis-rotation sequence of repeated extraction, so the gap cannot be uniquely attributed to CZ or measurement alone.
- Historical metadata remain archived, including the unused original unencoded r0 entry. That entry is excluded from the analysis.
- A record-availability timestamp is not a separately established logical-state preparation boundary. Server timestamps and full waveforms remain unavailable.
