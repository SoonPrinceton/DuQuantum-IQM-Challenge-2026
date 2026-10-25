"""Measurement conversion, intervals, optional ideal simulator and job collection."""

import json
import numpy as np
import stim

def extract(raw, reference):
    n = reference.num_measurements
    if set(raw) != {f'm_{i:03d}' for i in range(n)}:
        raise ValueError('Hardware measurement keys differ from the saved circuit')
    columns = [np.asarray(raw[f'm_{i:03d}'],dtype=bool) for i in range(n)]
    assert all(column.ndim==2 and column.shape[1]==1 for column in columns)
    measurements = np.concatenate(columns,axis=1)
    dets,obs = reference.compile_m2d_converter().convert(measurements=measurements,separate_observables=True)
    return measurements,dets,obs


def interval(k,n):
    z=1.959963984540054; p=k/n; den=1+z*z/n
    mid=(p+z*z/(2*n))/den; half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [max(0.,float(mid-half)),min(1.,float(mid+half))]


def ideal_raw(ops, shots):
    from qiskit import QuantumCircuit, transpile
    from qiskit_aer import AerSimulator
    names = sorted({q for op in ops for q in op['locus']})
    qmap = {q:i for i,q in enumerate(names)}
    n = sum(op['name']=='measure' for op in ops)
    qc = QuantumCircuit(len(names),n)
    for op in ops:
        qs = [qmap[q] for q in op['locus']]
        args = op['args']
        if op['name']=='prx':
            qc.rz(-args['phase'],qs[0]); qc.rx(args['angle'],qs[0]); qc.rz(args['phase'],qs[0])
        elif op['name']=='cz': qc.cz(*qs)
        elif op['name']=='barrier': qc.barrier(*qs)
        elif op['name']=='delay':
            for q in qs: qc.delay(args['duration'],q,unit='s')
        elif op['name']=='measure': qc.measure(qs[0],int(args['key'].split('_')[-1]))
    simulator = AerSimulator(method='matrix_product_state',max_parallel_threads=1)
    executable = transpile(qc,simulator,optimization_level=0)
    counts = simulator.run(executable,shots=shots,seed_simulator=20261024).result().get_counts()
    records = np.asarray([[int(c) for c in key.replace(' ','')[::-1]]
                          for key,count in counts.items() for _ in range(count)],dtype=bool)
    return {f'm_{i:03d}':records[:,i:i+1].astype(int).tolist() for i in range(n)}


def collect_job_results(job, expected_count, timeout_secs, status_path):
    """Native IQM CircuitJob.result() is nonblocking; wait explicitly first."""
    status = job.wait_for_completion(timeout_secs=timeout_secs)
    status_text = getattr(status, "value", str(status)).lower()
    data = getattr(job, "data", None)
    errors = [str(error) for error in getattr(data, "errors", [])]
    messages = [str(message) for message in getattr(data, "messages", [])]
    status_path.write_text(json.dumps({
        "job_id": str(job.job_id), "status": status_text,
        "errors": errors, "messages": messages,
    }, indent=2))
    if status_text != "completed":
        detail = "; ".join(errors) or "; ".join(messages) or "No server error details available"
        raise RuntimeError(
            f"Job {job.job_id} status={status_text}: {detail}. "
            "The job was not resubmitted or cancelled. Resume from this output directory."
        )
    raw = job.result()
    if raw is None or len(raw) != expected_count:
        received = "None" if raw is None else str(len(raw))
        raise RuntimeError(
            f"Completed job {job.job_id}: expected {expected_count} circuit results, got {received}. "
            "Job ID is saved; resume without resubmitting this batch."
        )
    return raw
