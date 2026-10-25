from surface_code import *
from stim_to_qiskit import stim_to_qiskit
from build_emerald_qubit_rotated import *
from extract_syndromes import *
from qiskit import transpile

import os
from dotenv import load_dotenv


load_dotenv()

TOKEN = os.getenv("IQM_TOKEN", None)
SERVER_URL = "https://resonance.iqm.tech"


def _strip_qrisp_model_annotations(qc):
    """Remove model-only instructions from a copy before hardware routing.

    Detector parities remain in DetectorExperiment's classical processing;
    noise annotations remain in its Stim model. Neither is a physical gate.
    Preserve all measurement instructions and classical-bit ordering.
    """
    from qrisp.jasp.primitives.parity_primitive import ParityOperation
    from qrisp.misc.stim_tools import StimNoiseGate

    hardware_qc = qc.copy()
    hardware_qc.data = [
        inst for inst in hardware_qc.data
        if not isinstance(inst.op, (ParityOperation, StimNoiseGate))
    ]
    return hardware_qc


# ─────────────────────────────────────────────────────────────────────────────
#  5.  HARDWARE EXECUTION ON IQM RESONANCE
#
#  Two equal tracks. Pick ONE. Both return the same syndromes dict
#  (at least "det_events" and "obs_flips"), so decode_hardware_results()
#  works for either.
#
#  NOTE: Emerald does not support mid-circuit reset -> use the no-reset circuit.
# ─────────────────────────────────────────────────────────────────────────────

def run_qiskit(
    distance: int    = 3,
    rounds: int      = 1,
    shots: int       = 5_000,
    token: str       = None,
    quantum_computer: str = "emerald",
    no_reset: bool   = True,
) -> dict:
    """
    Build in qiskit. Runs the surface code on IQM Emerald via IQM Resonance.

        Stim circuit → Qiskit circuit → transpile → run → counts
        → raw measurements → extract_syndromes()

    Parameters
    ----------
    distance         : code distance (use 3 for emerald)
    rounds           : measurement rounds
    shots            : number of circuit executions
    token            : IQM Resonance API token (falls back to TOKEN)
    quantum_computer : e.g. "emerald"
    no_reset         : build the no-reset circuit (required on Emerald)

    Returns
    -------
    syndromes : dict  output of extract_syndromes()
                "det_events" (shots, num_detectors), "obs_flips" (shots, num_observables), ...
    """
    from iqm.qiskit_iqm import IQMProvider

    # Build circuits, exchange for make_qiskit_circuit() if desired
    stim_circ  = make_stim_circuit(distance, rounds, no_reset=no_reset, noise=None)

    #CREATE A CIRCUIT TO BE EXECUTED ON HARDWARE
    qc = []
    #qc, stim_to_dense, meas_order = stim_to_qiskit(stim_circ)
    #qc = make_qiskit_circuit(args)


    # OPTIONAL: IF YOU WANT TO ASSIGN QUBIT MAP BY HAND
    # Map abstract Qiskit qubits to emerald physical qubits
    qubit_map  = build_emerald_qubit_map(stim_circ)

    # HIGHLY RECOMMEND PRINTING THIS, if you are using it.
    # print(qubit_map)
    # initial_layout: Qiskit needs a list where position i = physical qubit for logical i
    n_qubits   = qc.num_qubits
    layout     = [qubit_map.get(i, i) for i in range(n_qubits)]

    # Connect to IQM Resonance
    provider   = IQMProvider(
        SERVER_URL,
        quantum_computer=quantum_computer,
        token=token or TOKEN
    )
    backend    = provider.get_backend()

    # Transpile: Qiskit transpilerwill insert SWAPs for diagonal CX pairs automatically
    # think about whether this is necessary
    qc_t       = transpile(
    # arguments
    )
    print(f"Transpiled circuit depth: {qc_t.depth()}  "
          f"(ideal: {qc.depth()}, extra depth from SWAPs)")

    # Run
    job        = backend.run(qc_t, shots=shots)
    counts     = job.result().get_counts()

    # Convert counts dict → (shots, num_measurements) numpy array
    raw_meas   = counts_to_measurement_array(counts, stim_circ.num_measurements,
                                               shots)
    # m2d converter must come from the circuit that actually ran
    return extract_syndromes(raw_meas, stim_circ)



def run_qrisp(
    distance: int    = 3,
    rounds: int      = 1,
    shots: int       = 5_000,
    token: str       = None,
    quantum_computer: str = "emerald",
    no_reset: bool   = True,
) -> dict:
    """
    Build in qrisp. Runs the surface code on IQM Emerald via IQM Resonance.

        Stim circuit → Qrisp DetectorExperiment → IQMBackend → run
        → detection events + observable flips

    The DetectorExperiment stores the info the detectors and observables itself, so
    counts_to_measurement_array() and extract_syndromes() are not needed here.

    Parameters
    ----------
    distance         : code distance (use 3 for emerald)
    rounds           : measurement rounds
    shots            : number of circuit executions
    token            : IQM Resonance API token 
    quantum_computer : e.g. "emerald"
    no_reset         : build the no-reset circuit (required on Emerald)

    Returns
    -------
    syndromes : dict  same keys as extract_syndromes() from qiskit function
                "det_events" (shots, num_detectors), "obs_flips" (shots, num_observables)
    """
    import os
    from iqm.qrisp_iqm import IQMBackend
    from qrisp import CircuitPass
    from stim_to_qrisp import stim_to_detector_experiment

    os.environ["IQM_TOKEN"] = token or TOKEN

    # Build circuits, exchange for make_qrisp_experiment() if desired
    # stim_circ  = make_stim_circuit(distance, rounds, no_reset=no_reset, noise=None)

    # CREATE AN EXPERIMENT TO BE EXECUTED ON HARDWARE
    # experiment = []
    # experiment = stim_to_detector_experiment(stim_circ, no_reset=no_reset)
    experiment = make_qrisp_experiment(
        distance, rounds, no_reset=no_reset, p=0.003
    )

    # Connect to IQM Resonance
    # Layout: think about the pass_manager. iqm.qrisp_iqm.vf2pp_layout finds an
    # exact embedding on the coupling map (zero SWAPs) or fails.
    backend    = IQMBackend(
        server_url=SERVER_URL,
        device_instance=quantum_computer,
    # further arguments
    )

    # IQM's router assumes every instruction has at least one qubit.
    # Strip classical parity/model annotations before its default passes.
    backend.pm.insert_pass(0, CircuitPass(_strip_qrisp_model_annotations))

    # Run the experiment on the backend and collect detection events + observable flips
    # TO BE DONE

    model_circuit = experiment.to_stim()
    dem = model_circuit.detector_error_model(decompose_errors=True)

    det_events, obs_flips = experiment._sample_via_backend(
        shots=shots,
        backend=backend
    )

    det_events = np.array(det_events, dtype=bool)
    obs_flips  = np.array(obs_flips, dtype=bool)

    assert det_events.shape == (shots, model_circuit.num_detectors)
    assert obs_flips.shape == (shots, model_circuit.num_observables)

    return {"det_events": det_events, "obs_flips": obs_flips, "dem": dem, "shots":shots}



def decode_hardware_results(
    syndromes: dict
) -> tuple[float, float]:
    """
    Decodes the output of run_qiskit() or run_qrisp() via your desired decoder.

    NOTE: Other parameters can be important for this function, such as a noise model.
    If you decode with a detector error model (DEM): the DEM’s detectors must match the circuit that ran.


    Parameters
    ----------
    syndromes           : dict  output of run_qiskit() / run_qrisp()

    Returns
    -------
    ler : float   estimated logical error probability for this round count (fit over several round counts to get the logical error rate, see README)
    err : float   OPTIONAL 1-sigma statistical error
    """

    # USE e.g. PyMatching to decode

    matcher = pymatching.Matching.from_detector_error_model(
        syndromes["dem"]
    )

    pred = matcher.decode_batch(syndromes["det_events"])
    failed = np.any(pred != syndromes["obs_flips"], axis=1)

    ler = float(failed.mean())
    err = float(np.sqrt(ler * (1 - ler) / len(failed)))

    return ler, err

if __name__ == "__main__":
    syndromes = run_qrisp(distance=3, rounds=3, shots=1000, token=TOKEN, quantum_computer="emerald", no_reset=True)
    from pathlib import Path
    import numpy as np

    output = Path("diagnostics/results/surface_code/run_004_emerald")
    output.mkdir(parents=True, exist_ok=False)

    np.savez_compressed(
        output / "measurements.npz",
        det_events=syndromes["det_events"],
        obs_flips=syndromes["obs_flips"],
        shots=syndromes["shots"],
    )

    syndromes["dem"].to_file(str(output / "model.dem"))
    ler, err = decode_hardware_results(syndromes)
    print(f"Logical error probability: {ler:.4f} ± {err:.4f}")
