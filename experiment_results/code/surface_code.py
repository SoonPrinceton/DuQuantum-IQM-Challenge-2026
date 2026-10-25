"""Validated Stim generator and no-reset detector rewrite, extracted unchanged."""

from typing import Optional
import stim

DEFAULT_NOISE = dict(
    after_clifford_depolarization  = 0.003,
    after_reset_flip_probability   = 0.003,
    before_measure_flip_probability= 0.003,
    before_round_data_depolarization = 0.003,
)

def make_stim_circuit(
    distance: int,
    rounds: int,
    basis: str = "Z",
    no_reset: bool = True,
    noise: Optional[dict] = None,
) -> stim.Circuit:
    """Generate a rotated surface-code Stim reference (Z or X memory).

    Hardware uses the no-reset variant: remove ancilla resets while rewriting
    every detector/observable measurement expression. Only initial resets and
    their preparation errors remain. `noise=None` yields a noiseless reference.

    For the common unitary encoder and zero-round support used in the measured
    benchmark, see memory_benchmark.encoder() and memory_benchmark.model().
    """
    basis = basis.upper()

    if basis not in ("X", "Z"):
        raise ValueError("basis must be 'X' or 'Z'")
    if distance < 3:
        raise ValueError("distance must be at least 3")
    if rounds < 1:
        raise ValueError("rounds must be at least 1")

    noise_args = dict(noise or {})

    initial_reset_error = 0.0
    if no_reset:
        initial_reset_error = noise_args.pop(
            "after_reset_flip_probability", 0.0
        )
        noise_args["after_reset_flip_probability"] = 0.0

    circuit = stim.Circuit.generated(
        f"surface_code:rotated_memory_{basis.lower()}",
        distance=distance,
        rounds=rounds,
        **noise_args,
    )

    if no_reset:
        circuit = _remove_round_resets(
            circuit,
            initial_reset_error=initial_reset_error,
        )

    return circuit


def _remove_round_resets(
    circuit: stim.Circuit,
    initial_reset_error: float = 0.0,
) -> stim.Circuit:
    result = stim.Circuit()

    num_measured = 0

    last_measurement = {}

    measurement_terms = []

    for instruction in circuit.flattened():
        name = instruction.name
        targets = instruction.targets_copy()
        args = instruction.gate_args_copy()


        assert name not in ("MY", "MRX", "MRY", "MXX", "MZZ", "MPP")

        if name in ("M", "MR", "MX"):
            for target in targets:
                q = target.value

                terms = {num_measured}

                if name == "MR" and q in last_measurement:
                    terms.add(last_measurement[q])

                measurement_terms.append(terms)
                last_measurement[q] = num_measured
                num_measured += 1

            new_name = "M" if name == "MR" else name
            result.append(new_name, targets, args)

        elif name in ("DETECTOR", "OBSERVABLE_INCLUDE"):
            rewritten = set()

            for target in targets:
                original_index = num_measured + target.value

                rewritten ^= measurement_terms[original_index]

            new_targets = [
                stim.target_rec(index - num_measured)
                for index in sorted(rewritten)
            ]

            result.append(name, new_targets, args)

        else:
            result.append(instruction)

            if initial_reset_error:
                if name == "R":
                    result.append(
                        "X_ERROR", targets, initial_reset_error
                    )
                elif name == "RX":
                    result.append(
                        "Z_ERROR", targets, initial_reset_error
                    )

    return result
