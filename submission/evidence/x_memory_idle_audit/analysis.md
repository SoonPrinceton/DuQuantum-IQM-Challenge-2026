# Offline X-memory audit: 5 us versus 15 us

Scope: encoded logical |+_L>, **zero syndrome rounds** in both experiments.

**PASS:** preparation, full physical mapping, final readout and logical observable definitions are identical.
**PASS:** all 136 saved native instructions are identical except instruction 108: delay 5 us versus 15 us on the same 17 qubits.
Native counts: {'prx': 108, 'cz': 16, 'barrier': 2, 'delay': 1, 'measure': 9}. Logical X = X(QB1) X(QB4) X(QB9), evaluated through final X-basis data readout.
Both use calibration `23e12f6b-c23c-4736-8736-b0142989c0e6`.

| Saved pre-submission pulse timing | 5 us run | 15 us run |
|---|---:|---:|
| Initialization | 400 us | 400 us |
| Encoded state preparation | 0.752 us | 0.752 us |
| Memory idle | 5 us | 15 us |
| Final rotation/readout | 1.188 us | 1.188 us |
| Full schedule | 406.94 us | 416.94 us |

## Evidence limits

The source runs save native gate inputs and local compilation duration summaries, not the complete compiled pulse instructions, waveform samples, or downloaded server payloads.
Thus this audit establishes an idle-only difference at the native-input level and a consistent 10 us difference in the saved local schedule timings. It cannot prove waveform-by-waveform equivalence or the server-executed schedules.
An idle-only schedule change would shift all final-readout pulses by 10 us, so absolute pulse timestamps should not be compared without accounting for that shift.
Same calibration ID does not prove identical historical controller settings or prevent hardware drift between jobs.

## Decoder note

The noiseless reference, detector definitions and logical parity are identical. Noisy DEMs include storage-time-dependent idle noise, so their probabilities need not be identical; no pulse-equivalence claim follows from DEM equality.

See audit.json for exact hashes, the changed instruction, timing values and the four recorded job IDs.
