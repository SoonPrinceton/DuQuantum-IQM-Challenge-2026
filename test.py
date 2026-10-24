from surface_code import make_stim_circuit, DEFAULT_NOISE, get_circuit_info, make_qrisp_experiment
from raw_stim_d3_surface_code import logical_error_rate
import pymatching
import stim
import numpy as np
from stim_to_qrisp import stim_to_detector_experiment


def sanity_check_single_pauli(
    rounds: int = 3,
    no_reset: bool = True,
    shots: int = 16,
    verbose: bool = True,
) -> list[dict]:
    """Check each data-qubit X/Y/Z fault after the first ancilla readout.

    Uses a clean circuit's m2d reference, so a forced fault is not absorbed
    into the faulty circuit's reference sample. The decoder is built from
    DEFAULT_NOISE, independently of the forced fault. This checks a fixed
    time slice of a distance-3 memory circuit, not every circuit fault.
    """
    if shots < 1:
        raise ValueError("shots must be at least 1")

    results = []
    for basis in ("Z", "X"):
        clean = make_stim_circuit(3, rounds, basis, no_reset, noise=None)
        flat = list(clean.flattened())
        coords = clean.get_final_qubit_coordinates()
        data_qubits = sorted(
            q for q, xy in coords.items()
            if int(xy[0]) % 2 == 1 and int(xy[1]) % 2 == 1
        )
        ancillas = set(coords) - set(data_qubits)
        # Insert after the first complete ancilla measurement instruction.
        insertion = next(
            i + 1 for i, inst in enumerate(flat)
            if inst.name in ("M", "MR")
            and {t.value for t in inst.targets_copy()} == ancillas
        )
        prefix, suffix = stim.Circuit(), stim.Circuit()
        for inst in flat[:insertion]:
            prefix.append(inst)
        for inst in flat[insertion:]:
            suffix.append(inst)

        # Independently predict the raw logical parity from anticommutation
        # with the final logical readout support. No later gates rotate data
        # qubits in these generated memory circuits.
        measured_qubits = []
        logical_support = set()
        for inst in flat:
            if inst.name in ("M", "MR", "MX"):
                measured_qubits.extend(t.value for t in inst.targets_copy())
            elif inst.name == "OBSERVABLE_INCLUDE":
                assert inst.gate_args_copy() == [0.0]
                for t in inst.targets_copy():
                    logical_support ^= {measured_qubits[len(measured_qubits) + t.value]}
        assert logical_support and logical_support <= set(data_qubits)

        converter = clean.compile_m2d_converter()
        model = make_stim_circuit(3, rounds, basis, no_reset, DEFAULT_NOISE)
        matcher = pymatching.Matching.from_detector_error_model(
            model.detector_error_model(decompose_errors=True)
        )
        baseline_raw = clean.compile_sampler(seed=1234).sample(shots=shots)
        baseline_dets, baseline_obs = converter.convert(
            measurements=baseline_raw, separate_observables=True
        )
        assert not baseline_dets.any() and not baseline_obs.any()
        assert not matcher.decode_batch(baseline_dets).any()

        if verbose:
            print(f"\n{basis} memory: rounds={rounds}, no_reset={no_reset}")
            print("qubit  Pauli  fired detectors  actual observable  decoder prediction")
        for qubit in data_qubits:
            for pauli in ("X", "Y", "Z"):
                faulty = prefix.copy()
                # Probability 1 forces a physical error on every shot.
                faulty.append(f"{pauli}_ERROR", [qubit], 1.0)
                faulty += suffix
                raw = faulty.compile_sampler(seed=1234).sample(shots=shots)
                dets, obs = converter.convert(
                    measurements=raw, separate_observables=True
                )
                predicted = matcher.decode_batch(dets)
                expected_flip = qubit in logical_support and pauli != basis
                label = f"{basis} memory, {pauli} on qubit {qubit}"
                assert np.all(dets == dets[0]), label
                assert np.all(obs == expected_flip), f"Wrong observable: {label}"
                assert np.array_equal(predicted, obs), f"Decoder mismatch: {label}"
                fired = np.flatnonzero(dets[0]).tolist()
                row = dict(
                    basis=basis, qubit=qubit, pauli=pauli,
                    fired_detectors=fired, expected_flip=bool(expected_flip),
                    actual_flip=bool(obs[0, 0]),
                    predicted_flip=bool(predicted[0, 0]),
                )
                results.append(row)
                if verbose:
                    print(f"{qubit:5d}  {pauli:5s}  {str(fired):15s}  "
                          f"{int(obs[0, 0]):17d}  {int(predicted[0, 0]):18d}")
    if verbose:
        print(f"\nPassed {len(results)} single-Pauli cases ({shots} shots each).")
    return results

def sanity_check_no_noise(
        rounds: int = 3,
        no_reset: bool = True,
        shots: int = 1000,
        verbose: bool = True
):
    for basis in ("Z", "X"):
        circuit = make_stim_circuit(
            distance=3,
            rounds=rounds,
            basis=basis,
            no_reset=no_reset,
            noise=None,
        )

        dets, obs = circuit.compile_detector_sampler().sample(
            shots=shots,
            separate_observables=True,
        )

        assert not dets.any()
        assert not obs.any()

        circuit.detector_error_model()

        assert all(
            instruction.name != "MR"
            for instruction in circuit.flattened()
        )
    if verbose:
        print(f"Passed no-noise cases (rounds={rounds}, no_reset={no_reset}, shots={shots}).")

def sanity_check_stim2qrisp(rounds: int = 3):
    source = make_stim_circuit(
        distance=3,
        rounds=rounds,
        basis="Z",
        no_reset=True,
        noise=None,
    )
    experiment = stim_to_detector_experiment(source, no_reset=False)
    qc = experiment.to_qc()
    assert qc.num_qubits() == 17
    converted = experiment.to_stim()

    assert converted.num_measurements == source.num_measurements
    assert converted.num_detectors == source.num_detectors
    assert converted.num_observables == source.num_observables
        

if __name__ == "__main__":
    c = make_stim_circuit(3, 3, "Z", no_reset=True, noise=DEFAULT_NOISE)
    print(get_circuit_info(c))
    print(f"code distance (graphlike) = "
          f"{len(c.shortest_graphlike_error())}\n")
    sanity_check_stim2qrisp(rounds=3)
    sanity_check_no_noise(rounds=3, no_reset=True, shots=1000, verbose=True)
    sanity_check_single_pauli(rounds=3, no_reset=True, shots=16, verbose=True)