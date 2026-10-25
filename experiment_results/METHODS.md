# Implementation, reproduction and evidence details

A complete memory-experiment pipeline: Stim references → fixed physical layout
and native PRX/CZ circuits → IQM execution → measurement histories → detector
events → PyMatching minimum-weight perfect matching (MWPM).

This folder is self-contained. Its archived results are **byte-identical copies**
of the original measurements and metadata. Packaging performed only offline
checks: no new hardware jobs, experiments, parameter searches or simulated shots.

## Reproduce the analysis

Python 3.12 was used; exact observed package versions are recorded in
`environment_versions.json`. From this folder:

```bash
python -m pip install -r requirements.txt
python code/reproduce_analysis.py
```

The script works from other working directories too. No IQM account or network
access is needed once dependencies are installed. It checks archive hashes,
reconverts raw measurements, reproduces every saved MWPM prediction and error
count, and writes only derived outputs:

- `tables/summary.csv` and `tables/summary.json`: LER, raw error, shot counts,
  Wilson 95% intervals, calibration and job IDs.
- `tables/timing.csv`: requested and compiled memory windows, idle slots,
  initialization, encoding and readout costs, native gate counts.
- `plots/memory_15us.png`: logical X/Z error versus extraction rounds, including
  unencoded baselines and confidence intervals.
- `plots/zero_round_storage.png`: encoded zero-round memory at near-zero, 5 μs
  and 15 μs, with confidence intervals.
- `RESULTS.md` and `validation.json`: readable results and consistency checks.

To write derived outputs elsewhere, use
`python code/reproduce_analysis.py --output-dir /path/to/output`.
The script never modifies `data/`. `plots/memory_15us_original.png` is the original
archived plot; the other two plots are reproducible views of the same measurements.

## Experiment and timing definition

The rotated distance-3 code uses **9 data + 8 auxiliary qubits**, on one fixed
17-qubit Emerald patch. All experiments use calibration
`23e12f6b-c23c-4736-8736-b0142989c0e6`.

- Z memory prepares logical `|0_L⟩`; X memory prepares logical `|+_L⟩`.
- All encoded round counts use the same **unitary known-state preparation within
  each basis**, without ancilla measurement. Ancillas are uncomputed to `|0⟩`.
  This preparation is not fault tolerant or an arbitrary-input encoder.
- There are 0, 1, 2 or 3 **intermediate** syndrome-extraction rounds. No ancilla
  is reset between rounds. No physical correction gates are applied.
- Memory time starts **after initialization and encoded preparation**, and ends
  **before final basis rotation/readout**. Preparation/readout errors remain
  part of the measured end-to-end error probability even though their durations
  are outside the memory window.
- For target 15 μs, remaining idle time is split as evenly as possible over the
  `r+1` intervals before, between and after extraction. Actual extraction time is
  included. Saved local compilations give 15.000 μs for 0/3 rounds and the
  baselines, and 14.996 μs for 1/2 rounds: a reported 4 ns quantization shortfall.
- Unencoded `|0⟩` and `|+⟩` baselines use the same mapped data qubit, QB1,
  with the same memory-window definition. They receive no decoding correction.

With reset, raw measurements directly encode the syndrome. Without reset,
`m[k] = s[k] XOR m[k-1]`, so syndrome changes use `m[k] XOR m[k-2]`, with earlier
values initially zero. The known encoded initial state supplies deterministic
initial checks. Each encoded experiment retains all `8r` intermediate checks and
four final checks inferred from the final data readout, including **four checks
at zero rounds**. Each arm uses its own saved detector error model and complete
measurement history. The baseline decoder is ordinary MWPM, without correlations.

## Measured results

Each configuration has **2000 shots**, collected in two 1000-shot blocks.
The primary memory benchmark and its zero-round diagnostics comprise 14
configurations and 28,000 shots. The earlier four component controls add 8000
shots, giving **18 analyzed configurations and 36,000 measured shots** in total.
All calibration IDs, original job IDs and block ordering remain in `data/`.

| 15 μs target | Z-memory LER | X-memory LER |
|---|---:|---:|
| 0 extraction rounds | 25.10% | 48.50% |
| 1 extraction round | 32.30% | 45.45% |
| 2 extraction rounds | 36.50% | 49.90% |
| 3 extraction rounds | 37.40% | 48.60% |
| Unencoded single-qubit baseline | 0.65% | 29.55% |

| Encoded zero-round memory | Z-memory LER | X-memory LER |
|---|---:|---:|
| Near-zero storage, compiled 0 ns | 2.90% | 2.00% |
| 5.000 μs storage | 11.70% | 52.25% |
| 15.000 μs storage | 25.10% | 48.50% |

Confidence intervals and exact failure counts are in `RESULTS.md` and the tables.
These results do **not** demonstrate an advantage over the unencoded baseline
or improved protection from additional rounds under these conditions.
The low near-zero errors constrain the combined preparation/readout/decoder
floor, while degradation appears when storage is added. They do not isolate
preparation fidelity or establish a unique physical cause. X-memory is near
random logical parity at both 5 and 15 μs; these points do not support a precise
coherence-time fit. The duration groups were measured in separate jobs, so drift
remains possible despite the shared calibration ID.

## Timing and zero-SWAP evidence

`code/memory_benchmark.py` uses the existing IQM PulLA compiler, not gate-count
estimates. Globally synchronized prefix schedules locate preparation and memory
boundaries. Compiler-added terminal probes are removed from **local timing
prefixes only**, and their endpoints are inspected before idle-only schedule
cleanup. Full executable circuits retain the normal compiler and pulse
implementations. Initialization is 400 μs and is outside the memory window.

The saved X-memory 5/15 μs circuits have identical preparation, mapping and final
readout. Their 136 native instructions differ only at the one delay instruction.
Encoding is 0.752 μs, final rotation/readout is 1.188 μs, and full local schedule
lengths are 406.940/416.940 μs. The archived audit is in
`evidence/x_memory_idle_audit/`.

All later circuits match the direct no-routing PRX/CZ construction, and every CZ
edge is native in the saved architecture. This verifies no added routing gates
at the **native circuit level**, including decomposed-SWAP overhead relative to
that construction. `evidence/historical_zero_swap/` additionally preserves an
earlier Emerald job's verification against an independently unrouted reference
up to commuting gate reordering. That job used the earlier three-round Z-memory
circuit; it is not one of the 14 benchmark configurations.

**Evidence limit:** complete pulse waveforms, server schedules and historical
controller-setting snapshots were not archived. Local timing summaries cannot
prove waveform-by-waveform identity or server execution timing. The submission
does not claim those missing validations.

## Code and data layout

```text
code/                     Validated core, execution scripts and offline analysis
  surface_code.py         Stim generator and symbolic no-reset record rewrite
  geometry.py             Rotated-code geometry and check supports
  memory_benchmark.py     Common X/Z encoder, 15 μs construction, IQM/MWPM pipeline
  zero_round.py           Near-zero two-arm diagnostic, shared execution utility
  zero_round_5us.py        Same zero-round experiment with 5 μs storage
  utilities.py            Measurement conversion, intervals and job collection
  reproduce_analysis.py   Offline archive, mapping, decoding and timing checks
  component_controls.py   Earlier four-arm fixed-patch reconstruction and analysis
data/
  memory_15us/            Ten measured configurations
  memory_near_zero/       Two encoded zero-round configurations
  memory_5us/             Two encoded zero-round configurations
  source_patch/           Original layout-input plan and calibration architecture
  component_controls/     Four earlier full-schedule-matched Z-memory controls
evidence/                 Archived timing audit and historical zero-SWAP evidence
tables/                   Derived JSON/CSV
plots/                    Original and reproduced result plots
PROVENANCE.json           Source paths, copied-file hashes and code-copy changes
```

Per-arm data include raw measurement matrices, detector events, observed logical
flips and original decoder predictions in `block*.npz`; noiseless Stim references;
noisy Stim models; DEMs; and the native input instructions. Verbose duplicate raw
JSON, mock outputs, incomplete plans and exploratory parameter searches are
omitted. Original absolute paths inside archived metadata are preserved as
provenance only; the runnable code uses submission-local paths.

## Hardware implementation (optional; not needed for analysis)

The original execution and timing-compilation code is retained with portable
imports and paths. Hardware dependencies are listed separately:

```bash
python -m pip install -r requirements-hardware.txt
python code/memory_benchmark.py --help
python code/zero_round.py --help
python code/zero_round_5us.py --help
```

For a future authorized rerun, provide `IQM_TOKEN` through the environment or a
local `.env` (no credentials are included). Each script's `--prepare` compiles
and validates a **new** plan under `output/`; it does not submit jobs. Only an
explicit `--hardware --resume-dir NEW_PREPARED_DIRECTORY` submits. Job IDs and
raw results are checkpointed to avoid duplicate submission on resume. Do not
use the archival `data/` folders as new run output directories.

Optional `--self-test` checks include noiseless Stim/native equivalence and
single-Pauli recovery; packaging did not run new shot-sampling experiments.
The default effective noise model is not calibrated T1/T2: gate/preparation/
readout probabilities are 0.003, with a declared effective depolarization model
for timed idle. A per-round probability must not be interpreted as a per-μs
rate. No pulse implementation is modified.

Source acknowledgements and exact extraction details are in `PROVENANCE.json`.
The original project and `diagnostics/` remain unchanged.

## Earlier component comparison

`COMPONENT_CONTROLS.md` and `plots/component_controls.png` compare 1 round + idle,
1 round + 48 extra CZ echo gates, 1 round + 16 extra ancilla measurements and
3 rounds. These are the earlier Z-memory circuits with their original preparation
convention and primary DEMs; they are not redefined as the later common unitary
encoder or as 15 μs memory-window experiments.

Their full local schedule durations are 410.872, 410.876, 410.872 and 410.872 μs;
padding is respectively 8.240, 7.904, 6.168 and 0 μs. Isolated preparation/memory/
readout boundaries were not archived, so their CSV fields remain blank. All
physical CZ edges match the saved calibration, and input sequences are reproduced
from the original fixed-patch source plus the specified interventions.

The one-round primary DEMs are identical. Extra CZs are ideal identity pairs;
extra Z measurements read already-measured ancillas without H or data coupling.
Those 16 auxiliary bits are archived in their original named-key JSON, checked
against the main measurement matrices and excluded from the decoder by design.

All four original failure counts, per-shot predictions and Wilson intervals are
reproduced. Derived contrasts use independent shots, Newcombe/Wilson difference
intervals and Fisher tests with Holm correction across the **three displayed
contrasts**. The original five-arm metadata/statistics remain unchanged; their
historical contrast family also included the omitted unencoded zero-round arm.
Non-significant echo/readout contrasts are not equivalence demonstrations.
