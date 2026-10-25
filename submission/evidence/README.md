# Evidence scope

- `historical_zero_swap/`: byte-identical archived gate-level payload,
  calibration architecture, execution parameters and the prior verification
  for Emerald job `01a16fc6-72f2-7186-9112-3ea5083a060f`. The stored verification
  reports equivalence to a reference compiled without layout/routing, up to
  commuting reordering and global phase. It does not claim historical pulse
  compiler logs or waveforms were reconstructed.
- `x_memory_idle_audit/`: byte-identical offline audit of the measured 5/15 μs
  zero-round logical X circuits. It records exact preparation/readout hashes,
  mapping, the changed delay instruction, local timing summaries and job IDs.
  Its original source paths are provenance, not runtime dependencies.
- The primary datasets' `pre_submission_timing.json` files contain their
  archived local pulse-compilation duration records, including initialization,
  preparation, storage and readout decomposition.

`code/reproduce_analysis.py` independently checks every primary circuit against
the submission's direct no-routing construction and saved native CZ topology,
and verifies the saved X-memory 5/15 μs idle-only input difference.

Full pulse instructions/waveforms, server schedules and historical controller
settings are missing. They have not been fabricated or inferred from gate counts.
