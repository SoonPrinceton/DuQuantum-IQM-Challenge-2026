# Distance-3 surface-code memory on IQM Emerald

An end-to-end pipeline for rotated surface-code memory: encoded logical X/Z
states, no-reset syndrome extraction, a fixed zero-routing-overhead patch,
IQM native hardware execution and PyMatching MWPM decoding.

This self-contained submission contains **18 measured configurations, 2000 shots
each (36,000 shots)**. Archived measurements and metadata are byte-identical to
the originals. No new experiments or hardware jobs were run to prepare it.

## Reproduce offline

Python 3.12 was used. From this folder:

```bash
python -m pip install -r requirements.txt
python code/reproduce_analysis.py
```

No IQM account or network is needed after dependencies are installed. This
reconverts the archived raw measurements, reproduces every original MWPM
prediction, checks counts/intervals/mapping/timing, and generates:

- `tables/summary.csv`, `tables/summary.json`, `tables/timing.csv`
- `plots/memory_15us.png`, `plots/zero_round_storage.png`
- `tables/component_controls.json`, `plots/component_controls.png`
- `RESULTS.md`, `validation.json`

The script never changes `data/`. Use `--output-dir PATH` for alternative derived
output. `plots/memory_15us_original.png` preserves the original archived plot.

## Experimental results

The 15 μs target is the window **after encoded preparation and before final
basis rotation/readout**, including extraction and idle. Idle is distributed
over `r+1` intervals. Actual compiled windows are 15.000 μs for 0/3 rounds and
the baselines, and 14.996 μs for 1/2 rounds (4 ns quantization shortfall).

| Intermediate syndrome rounds | Z-memory LER | X-memory LER |
|---|---:|---:|
| 0 | 25.10% | 48.50% |
| 1 | 32.30% | 45.45% |
| 2 | 36.50% | 49.90% |
| 3 | 37.40% | 48.60% |
| Unencoded single-qubit baseline | 0.65% | 29.55% |

The zero-round diagnostics preserve the **same complete encoder and readout**:

| Storage | Z-memory LER | X-memory LER |
|---|---:|---:|
| Near zero, compiled 0 ns | 2.90% | 2.00% |
| 5.000 μs | 11.70% | 52.25% |
| 15.000 μs | 25.10% | 48.50% |

Exact counts and Wilson 95% intervals are in `RESULTS.md` and the tables.
These measurements **do not demonstrate error-correction advantage** over the
unencoded baseline. Near-zero errors constrain the combined preparation/readout/
decoder floor; they do not measure preparation fidelity alone. Added storage
substantially degrades memory, especially X-memory. Separate jobs and possible
drift prevent a unique physical diagnosis or a precise coherence-time fit.

### Component controls: rerun with matched first-record windows

| Control | Failures / shots | Decoded LER | Raw parity error |
|---|---:|---:|---:|
| 1 round + idle | 279/2000 | 13.95% | 32.30% |
| 1 round + CZ echo | 285/2000 | 14.25% | 30.85% |
| 1 round + repeated ancilla measurement | 250/2000 | 12.50% | 29.50% |
| 3 rounds (historical) | 651/2000 | 32.55% | 41.70% |

The three one-round arms were remeasured with shorter padding, matching the historical three-round **first complete syndrome record → final data-readout barrier** interval to approximately 5.348 μs.
The three-round measurement data are reused, not remeasured. These are different acquisitions; shared calibration does not eliminate drift.
**This is record-window matching, not a verified common independent encoded-state preparation boundary.** The three-round circuit is pipelined.
The previous full-schedule-matched data are preserved under `evidence/previous_component_controls/`. Full original metadata for both acquisitions are in `data/component_controls/acquisition_sets/`.
See [COMPONENT_CONTROLS.md](COMPONENT_CONTROLS.md) for confidence intervals, timing and contrasts.

## Implementation and evidence

- Distance 3: 9 data + 8 auxiliary qubits. All encoded round counts use a common
  unitary known-state encoder within each basis. It is not fault tolerant.
- No intermediate resets or physical correction gates. No-reset detection
  events use `m[k] XOR m[k-2]`; all available syndrome history and four final
  data-derived checks are decoded with each arm's saved DEM.
- Calibration: `23e12f6b-c23c-4736-8736-b0142989c0e6`. Every configuration has
  two 1000-shot blocks; original job IDs, ordering and parameters are retained.
- Native sequences match the direct fixed-layout no-routing construction and
  every CZ link is present in the archived architecture. Historical independent
  zero-SWAP verification is preserved in `evidence/historical_zero_swap/`.
- Timings come from the existing **IQM pulse compiler**, not gate counts.
  Archived preparation/memory/readout boundaries are checked offline.
- **Missing evidence is explicit:** complete waveforms, server schedules and
  historical controller settings were not saved. Gate-level equivalence and
  local duration records do not prove server pulse-by-pulse equivalence.

## Folder guide

| Path | Contents |
|---|---|
| `code/` | Stim/no-reset core, geometry, benchmark and diagnostic execution, offline analysis |
| `data/memory_15us/` | Ten configurations: X/Z, 0–3 rounds, single-qubit baselines |
| `data/memory_near_zero/`, `data/memory_5us/` | Two encoded zero-round configurations each |
| `data/component_controls/` | Four requested earlier Z-memory component controls |
| `data/source_patch/` | Original layout-input plan and calibration architecture |
| `evidence/` | Historical zero-SWAP and 5/15 μs X-memory audit records |
| `tables/`, `plots/` | Derived tables and figures |
| `PROVENANCE.json` | Original source paths, copied-file hashes and code-copy changes |

Per-arm archives include raw measurement matrices, detector events, logical
flips and predictions (`block*.npz`), Stim circuits, DEMs and native instructions.
Duplicate raw JSON, debug/mock outputs and exploratory searches are omitted.
Original absolute metadata paths remain as provenance; runtime paths are local.
Component acquisition metadata are preserved separately per source. The old unused unencoded zero-round entry remains only in historical provenance and is excluded from analyses.

For implementation details and optional future hardware execution, see
**[METHODS.md](METHODS.md)**. Hardware dependencies are separate in
`requirements-hardware.txt`; only explicit `--hardware` submits jobs. Credentials
are not included. The original project and `diagnostics/` remain unchanged.
