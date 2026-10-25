"""15 us Emerald quantum memory AFTER initialization and state preparation.

python code/memory_benchmark.py --self-test
python code/memory_benchmark.py --prepare
python code/memory_benchmark.py --hardware --resume-dir PREPARED_DIR
python code/reproduce_analysis.py

Preparation reads calibration but never submits. The timing window starts at
the global barrier after encoded state preparation and ends at the global
barrier BEFORE final basis rotation/readout. Pulse-scheduled prefixes determine
both boundaries; no guessed reset time or gate-count estimate is subtracted.
No pulse implementation or controller settings are overridden.
"""

import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
from dataclasses import replace
import json
from itertools import combinations
import math
import os
from pathlib import Path
import sys
from uuid import UUID

import numpy as np
import pymatching
import stim
from surface_code import make_stim_circuit, DEFAULT_NOISE
from geometry import description
from utilities import extract, interval, ideal_raw
from utilities import collect_job_results

HERE = Path(__file__).resolve().parent
TARGET_NS = 15_000.
TOLERANCE_NS = 4.  # Controller sample-grid tolerance; never accept an overrun.
SOURCE = HERE.parent / 'data/source_patch'


def logical_support(basis):
    source = make_stim_circuit(3, 1, basis=basis, no_reset=True)
    record = [t.value for i in source.flattened() if i.name in ('M', 'MX')
              for t in i.targets_copy()]
    return [record[len(record) + t.value] for i in source.flattened()
            if i.name == 'OBSERVABLE_INCLUDE' for t in i.targets_copy()]


def encoder(basis):
    """Known-state, unitary encoder; all eight bus ancillas end in |0>.

    For X memory prepare the dual-code orbit, then H on every data qubit.
    This is not the transversal H of the original code without a permutation.
    """
    d = description()
    chosen = set(d['x_ancillas']) if basis == 'Z' else set(d['ancillas']) - set(d['x_ancillas'])
    frequency = Counter(q for a in chosen for q in d['checks'][a])
    c = stim.Circuit(); c.append('R', sorted(d['coords']))
    pivots = {a: min(q for q in d['checks'][a] if frequency[q] == 1) for a in sorted(chosen)}
    c.append('H', list(pivots.values()))
    for a, pivot in pivots.items():
        c.append('CX', [pivot, a])
        for q in sorted(d['checks'][a] - {pivot}): c.append('CX', [a, q])
        c.append('CX', [pivot, a])
    if basis == 'X': c.append('H', d['data'])
    return c


def round_body():
    """Reuse the project's extraction order; measured ancillas are NOT reset."""
    c = stim.Circuit()
    for inst in make_stim_circuit(3, 1, basis='Z', no_reset=True).flattened():
        if inst.name in ('H', 'CX', 'TICK'): c.append(inst)
        elif inst.name == 'M':
            c.append(inst)
            break  # First eight ancilla measurements, before final data readout.
    assert c.num_measurements == 8
    return c


def check_encoded_state(s, basis):
    d = description()
    def check(support, axis):
        p = stim.PauliString(max(d['coords']) + 1)
        for q in support: p[q] = axis
        assert s.peek_observable_expectation(p) == 1, (basis, support, axis)
    for a, support in d['checks'].items(): check(support, 'X' if a in d['x_ancillas'] else 'Z')
    check(logical_support(basis), basis)
    for a in d['ancillas']: check([a], 'Z')


def preparation_check(basis):
    s = stim.TableauSimulator(); s.do_circuit(encoder(basis))
    check_encoded_state(s, basis)


def native_preparation_check(arm, layout):
    """Check all stabilizers BEFORE storage, including native phase conventions."""
    inverse = {v: q for q, v in layout.items()}; s = stim.TableauSimulator()
    cache = {}
    for inst in arm['prep_ops']:
        qs = [inverse[q] for q in inst['locus']]
        if inst['name'] == 'cz': s.cz(*qs)
        else:
            angle, phase = inst['args']['angle'], inst['args']['phase']
            key = angle, phase
            if key not in cache:
                axis = np.array([[0., np.exp(-1j*phase)], [np.exp(1j*phase), 0.]])
                u = math.cos(angle/2)*np.eye(2) - 1j*math.sin(angle/2)*axis
                cache[key] = stim.Tableau.from_unitary_matrix(u, endian='little')
            s.do_tableau(cache[key], qs)
    check_encoded_state(s, arm['basis'])


def with_gate_noise(circuit, enabled):
    result = stim.Circuit()
    for inst in circuit:
        result.append(inst)
        if not enabled: continue
        qs = [t.value for t in inst.targets_copy()]
        if inst.name == 'R': result.append('X_ERROR', qs, DEFAULT_NOISE['after_reset_flip_probability'])
        elif inst.name == 'H': result.append('DEPOLARIZE1', qs, DEFAULT_NOISE['after_clifford_depolarization'])
        elif inst.name == 'CX': result.append('DEPOLARIZE2', qs, DEFAULT_NOISE['after_clifford_depolarization'])
    return result


def model(arm, noise=False, idle_error_per_us=.003, fault=None):
    """Full history, including known initial checks and final data-derived checks.

    Idle noise is an effective Pauli model, not calibrated T1/T2. It acts on
    live ancillas too, and is inserted at EACH of the r+1 hardware idle slots.
    fault=(qubit, Pauli) injects a terminal data fault for validation only.
    """
    d = description(); basis = arm['basis']; r = arm['rounds']; single = arm['kind'] == 'single'
    c = stim.Circuit()
    for q, xy in d['coords'].items():
        if not single or q == arm['data_qubit']: c.append('QUBIT_COORDS', [q], xy)
    if single:
        prep = stim.Circuit(); prep.append('R', [arm['data_qubit']])
        if basis == 'X': prep.append('H', [arm['data_qubit']])
    else: prep = encoder(basis)
    c += with_gate_noise(prep, noise)
    anc_order = [t.value for i in round_body() if i.name == 'M' for t in i.targets_copy()]
    data = [arm['data_qubit']] if single else d['data']
    active = data if single else sorted(d['coords'])
    for slot, duration in enumerate(arm['idle_ns']):
        c.append('TICK')
        if noise and duration > 0 and idle_error_per_us > 0:
            p = .75 * (1 - (1 - 4 * idle_error_per_us / 3) ** (duration / 1000.))
            c.append('DEPOLARIZE1', active, p)
        if slot == r: break
        if noise: c.append('DEPOLARIZE1', data, DEFAULT_NOISE['before_round_data_depolarization'])
        for inst in with_gate_noise(round_body(), noise):
            if inst.name == 'M' and noise:
                c.append('M', inst.targets_copy(), DEFAULT_NOISE['before_measure_flip_probability'])
            else: c.append(inst)
        # Raw m[k] = syndrome[k] XOR m[k-1]. Thus changes of syndrome
        # are m[k] XOR m[k-2], with both earlier values initially zero.
        n = 8 * (slot + 1)
        for j, a in enumerate(anc_order):
            terms = [stim.target_rec(8 * slot + j - n)]
            if slot >= 2: terms.append(stim.target_rec(8 * (slot - 2) + j - n))
            c.append('DETECTOR', terms, [*d['coords'][a][:2], slot])
    if fault: c.append(fault[1], [fault[0]])
    c.append('MX' if basis == 'X' else 'M', data,
             DEFAULT_NOISE['before_measure_flip_probability'] if noise else 0.)
    n = 8 * r + len(data); positions = {q: 8 * r + j for j, q in enumerate(data)}
    if not single:
        for a in d['ancillas']:
            if (a in d['x_ancillas']) != (basis == 'X'): continue
            terms = [stim.target_rec(positions[q] - n) for q in sorted(d['checks'][a])]
            if r:
                j = anc_order.index(a)
                terms.append(stim.target_rec(8 * (r - 1) + j - n))
                if r >= 2: terms.append(stim.target_rec(8 * (r - 2) + j - n))
            c.append('DETECTOR', terms, [*d['coords'][a][:2], r])
    logical = data if single else logical_support(basis)
    c.append('OBSERVABLE_INCLUDE', [stim.target_rec(positions[q] - n) for q in logical], 0)
    assert c.num_detectors == (0 if single else 8 * r + 4)
    return c


def operation(name, locus, **args):
    return dict(name=name, locus=list(locus), args=args)


def lower_unitary(circuit, layout):
    """Exact Clifford -> native gates on the fixed patch, with zero routing.

    H = RY(-pi/2) RX(pi), up to global phase. CX = H(target) CZ H(target).
    Unlike a terminal-Z optimization, this preserves the state at idle slots.
    PulLA's existing compiler supplies all physical pulse implementations.
    """
    out = []
    def h(q):
        out.extend([operation('prx', [layout[q]], angle=float(math.pi), phase=0.),
                    operation('prx', [layout[q]], angle=float(-math.pi / 2), phase=float(math.pi / 2))])
    for inst in circuit:
        qs = [t.value for t in inst.targets_copy()]
        if inst.name == 'R': continue  # Standard compiler supplies initial ground-state preparation.
        if inst.name == 'H':
            for q in qs: h(q)
        elif inst.name == 'CX':
            for a, b in zip(qs[::2], qs[1::2]):
                h(b); out.append(operation('cz', [layout[a], layout[b]])); h(b)
        elif inst.name != 'TICK': raise ValueError(inst.name)
    return out


def memory_operations(arm):
    out = []; locus = arm['active_locus']
    for slot, duration in enumerate(arm['idle_ns']):
        out.append(operation('barrier', locus))
        if duration > 0: out.append(operation('delay', locus, duration=float(duration * 1e-9)))
        out.append(operation('barrier', locus))
        if slot < arm['rounds']:
            out += deepcopy(arm['round_ops'])
            out.extend(operation('measure', [q], key=f'm_{8 * slot + j:03d}')
                       for j, q in enumerate(arm['ancilla_locus']))
    return out


def instructions(arm):
    return deepcopy(arm['prep_ops']) + memory_operations(arm) + deepcopy(arm['readout_ops'])


def make_arms(layout, baseline_qubit=None):
    d = description(); body = round_body()
    anc = [t.value for i in body if i.name == 'M' for t in i.targets_copy()]
    unitary = stim.Circuit()
    for i in body:
        if i.name != 'M': unitary.append(i)
    baseline_qubit = baseline_qubit or layout[d['data'][0]]
    inverse = {v: k for k, v in layout.items()}
    if baseline_qubit not in {layout[q] for q in d['data']}:
        raise ValueError('Baseline must use a data qubit on the fixed patch')
    arms = []
    for basis in ('Z', 'X'):
        preparation_check(basis)
        for kind, r in [('encoded', r) for r in range(4)] + [('single', 0)]:
            data = d['data'] if kind == 'encoded' else [inverse[baseline_qubit]]
            prep = encoder(basis) if kind == 'encoded' else stim.Circuit(f'R {data[0]}')
            if kind == 'single' and basis == 'X': prep.append('H', data)
            readout = stim.Circuit()
            if basis == 'X': readout.append('H', data)
            name = f'{basis}_r{r}' if kind == 'encoded' else f'{basis}_single'
            arm = dict(name=name, basis=basis, kind=kind, rounds=r, data_qubit=data[0],
                       active_locus=[layout[q] for q in sorted(d['coords'])] if kind == 'encoded' else [baseline_qubit],
                       prep_ops=lower_unitary(prep, layout), round_ops=lower_unitary(unitary, layout),
                       ancilla_locus=[layout[q] for q in anc], idle_ns=[0.] * (r + 1),
                       readout_ops=lower_unitary(readout, layout))
            arm['readout_ops'] += [operation('measure', [layout[q]], key=f'm_{8*r+j:03d}') for j, q in enumerate(data)]
            if kind == 'encoded' and r == 0: native_preparation_check(arm, layout)
            arms.append(arm)
    return arms


def native_circuits(arms):
    from iqm.pulse import Circuit, CircuitOperation
    return [Circuit(name=a['name'], instructions=tuple(CircuitOperation(
        name=i['name'], locus=tuple(i['locus']), args=i['args']) for i in instructions(a))) for a in arms]


def validate(arms, native=True):
    """Noiseless Stim, terminal single-Pauli recovery, and native-vs-Stim checks."""
    d = description()
    for a in arms:
        ref = model(a); dets, obs = ref.compile_detector_sampler(seed=17).sample(128, separate_observables=True)
        assert not dets.any() and not obs.any(), a['name']
        dem = model(a, True).detector_error_model(decompose_errors=True)
        matcher = pymatching.Matching.from_detector_error_model(dem) if a['kind'] == 'encoded' else None
        for q in d['data'] if matcher else [a['data_qubit']]:
            fault = model(a, fault=(q, 'X' if a['basis'] == 'Z' else 'Z'))
            # Convert faulted measurements with the CLEAN reference. Sampling a
            # faulted circuit's own detectors would subtract its fault reference.
            raw = fault.compile_sampler(seed=23).sample(4)
            ds, os_ = ref.compile_m2d_converter().convert(measurements=raw, separate_observables=True)
            if matcher: assert np.array_equal(matcher.decode_batch(ds), os_), (a['name'], q)
            else: assert os_.all()
        if native:
            _, ds, os_ = extract(ideal_raw(instructions(a), 32), ref)
            assert not ds.any() and not os_.any(), ('native', a['name'])
    print('PASS: encoded X/Z preparation, noiseless Stim, complete no-reset history, single-Pauli recovery, native readout.', flush=True)


def compile_durations(compiler, arms):
    from iqm.pulse import Circuit, CircuitOperation
    # Compile prefixes alongside full circuits, with identical global boundary
    # barriers, calibration, settings and initialization context. The global
    # barriers serialize prep/memory/readout on the complete active patch.
    # Prefixes are compiler probes ONLY: they are never submitted to hardware.
    circuits = []
    for a in arms:
        barrier = [operation('barrier', a['active_locus'])]
        prefixes = [barrier, a['prep_ops'] + barrier,
                    a['prep_ops'] + memory_operations(a), instructions(a)]
        for label, ops in zip(('initialization', 'prepared', 'memory_end', 'full'), prefixes):
            circuits.append(Circuit(name=f'{a["name"]}_{label}', instructions=tuple(
                CircuitOperation(name=i['name'], locus=tuple(i['locus']), args=i['args']) for i in ops)))
    settings = compiler.get_settings(circuits=circuits); settings.set_shots(1)
    context = compiler.compiler_context(None, settings); context['timebox_input'] = False
    circuits, context = compiler.run_stages(circuits, context, stages=compiler.circuit_stages)
    from iqm.cpc.compiler.standard_stages import ONLY_PROBE_KEY, ADDITIONAL_SUBSCRIPTION_KEY
    for j, circuit in enumerate(circuits):
        if j % 4 != 3:
            # subscribe_and_probe adds terminal readout even to a measurement-
            # free prefix. Remove ONLY those automatic additions from compiler
            # probes. Full executable circuits retain the default readout and
            # pulse stages unchanged. Explicit syndrome measurements remain.
            circuits[j] = replace(circuit, instructions=tuple(i for i in circuit.instructions
                if i.args.get('key') not in (ONLY_PROBE_KEY, ADDITIONAL_SUBSCRIPTION_KEY)))
    stages = list(compiler.pulse_stages)
    split = next(i for i, stage in enumerate(stages) if stage.name == 'schedule_stage')
    schedules, context = compiler.run_stages(circuits, context, stages=stages[:split])
    spot = context['spot_contexts'][0]; builder = spot['builder']; metrics = spot.get('circuit_metrics')
    # Inspect scheduled prefix endpoints BEFORE cleanup: an initialization-only
    # or idle-only probe has exclusively Waits, which cleanup correctly removes
    # from executable programs but which must be retained for boundary timing.
    durations = [max(builder.channels[ch].duration_to_seconds(seg.duration) for ch, seg in schedule.items()) * 1e9
                 for schedule in schedules]
    finished, context = compiler.run_stages(schedules, context, stages=stages[split:])
    full_durations = [max(builder.channels[ch].duration_to_seconds(seg.duration) for ch, seg in finished[4*j+3].items()) * 1e9
                      for j in range(len(arms))]
    rows = []
    for j, a in enumerate(arms):
        init, prepared, memory_end, full = durations[4*j:4*j+4]
        assert abs(full-full_durations[j]) <= TOLERANCE_NS, 'Cleanup changed full-schedule boundary; cannot use prefix timing'
        full = full_durations[j]
        assert init <= prepared + 1e-6 and prepared <= memory_end + 1e-6 and memory_end <= full + 1e-6, 'Non-monotone prefix schedules'
        assert full - memory_end > 0, 'Final readout must be outside the memory boundary'
        include_readout = a.get('include_final_readout', False)
        window = (full if include_readout else memory_end) - prepared
        row = dict(name=a['name'], total_schedule_ns=float(full), idle_ns=a['idle_ns'].copy(),
                   initialization_ns=float(init), state_preparation_ns=float(prepared-init),
                   storage_before_readout_ns=float(memory_end-prepared), final_readout_ns=float(full-memory_end),
                   benchmark_duration_ns=float(window), includes_final_readout=include_readout,
                   prepared_boundary_ns=float(prepared), readout_boundary_ns=float(memory_end))
        if metrics:
            cz = int(sum(sum(c.values()) for c in metrics[4*j+3].gate_loci.get('cz', {}).values()))
            row['compiled_cz_count'] = cz
            assert cz == sum(i['name'] == 'cz' for i in instructions(a)), 'Compiler changed CZ count'
        rows.append(row)
    return rows


def fit_timing(compiler, arms, compile_fn=compile_durations):
    """Fit the post-preparation pulse window, accounting for nonzero round time."""
    initial = compile_fn(compiler, arms); history = [initial]
    eligible = []
    for a, row in zip(arms, initial):
        a['unpadded_schedule_ns'] = row['total_schedule_ns']
        a['unpadded_duration_ns'] = row['benchmark_duration_ns']
        a['fits'] = row['benchmark_duration_ns'] <= TARGET_NS
        if a['fits']:
            # Rounds and all other nonzero costs were already pulse-compiled.
            # Split remaining time on a 4 ns grid, differing by at most one tick.
            ticks = int(math.floor((TARGET_NS - row['benchmark_duration_ns']) / 4 + 1e-8))
            n = a['rounds'] + 1; q, rem = divmod(ticks, n)
            a['idle_ns'] = [float(4 * (q + (i < rem))) for i in range(n)]
            eligible.append(a)
        else: a['timing_failure'] = f'Unpadded post-preparation window {row["benchmark_duration_ns"]:.1f} ns exceeds 15000 ns'
    # Adding a delay can change scheduler alignment. Recompile every iteration.
    for _ in range(8):
        if not eligible: break
        rows = compile_fn(compiler, eligible); history.append(rows); pending = []
        for a, row in zip(eligible, rows):
            a['verified_schedule_ns'] = row['total_schedule_ns']
            a['verified_duration_ns'] = row['benchmark_duration_ns']
            a['timing_details'] = row
            delta = TARGET_NS - row['benchmark_duration_ns']
            if -1e-6 <= delta <= TOLERANCE_NS: continue
            total = sum(a['idle_ns']) + delta
            if total < 0:
                a['fits'] = False; a['timing_failure'] = 'Scheduler alignment cannot fit'; continue
            ticks = max(0, int(math.floor(total / 4 + 1e-8))); n = a['rounds'] + 1
            q, rem = divmod(ticks, n); a['idle_ns'] = [float(4 * (q + (i < rem))) for i in range(n)]
            pending.append(a)
        if not pending: break
        eligible = pending
    # Equal slots can quantize differently depending on where the longer slot
    # sits relative to extraction. Try the few balanced orderings near target
    # before accepting a one-tick shortfall (at most four slots per arm).
    probes = []; owners = []
    for a in arms:
        delta = TARGET_NS - a.get('verified_duration_ns', -1.)
        if not a['fits'] or not 1e-6 < delta <= TOLERANCE_NS: continue
        n = a['rounds'] + 1
        for extra_ticks in (0, 1, 2):
            ticks = round(sum(a['idle_ns']) / 4) + extra_ticks; q, rem = divmod(ticks, n)
            for longer in combinations(range(n), rem):
                trial = deepcopy(a); trial['idle_ns'] = [float(4*(q+(i in longer))) for i in range(n)]
                if trial['idle_ns'] == a['idle_ns']: continue
                probes.append(trial); owners.append(a)
    if probes:
        rows = compile_fn(compiler, probes); history.append(rows)
        for trial, a, row in zip(probes, owners, rows):
            delta = TARGET_NS - row['benchmark_duration_ns']
            if -1e-6 <= delta < TARGET_NS - a['verified_duration_ns'] - 1e-6:
                a['idle_ns'] = trial['idle_ns']; a['timing_details'] = row
                a['verified_schedule_ns'] = row['total_schedule_ns']; a['verified_duration_ns'] = row['benchmark_duration_ns']
    for a in arms:
        ns = a.get('verified_duration_ns', a['unpadded_duration_ns'])
        a['fits'] = bool(a['fits'] and -1e-6 <= TARGET_NS - ns <= TOLERANCE_NS)
        if not a['fits']: a.setdefault('timing_failure', 'Timing fit did not converge')
    increments = {}
    for basis in ('X', 'Z'):
        group = sorted((a for a in arms if a['kind'] == 'encoded' and a['basis'] == basis), key=lambda a: a['rounds'])
        increments[basis] = [dict(rounds=a['rounds'], added_schedule_ns=a['unpadded_duration_ns']-b['unpadded_duration_ns'])
                             for b, a in zip(group, group[1:])]
    return dict(target_ns=TARGET_NS, tolerance_ns=TOLERANCE_NS, includes_initialization=False,
                includes_state_preparation=False, includes_final_readout=arms[0].get('include_final_readout', False),
                prefix_terminal_probes_removed=True,
                method='Compiled global-barrier prefixes in the same batch as full schedules', server_timestamp_verified=False,
                all_fit=all(a['fits'] for a in arms), incremental_round_cost_ns=increments, iterations=history)


def timing_self_test(arms):
    """Exercise fitting and rejection without a network or hardware job."""
    def fake_compile(_, group):
        return [dict(name=a['name'], total_schedule_ns=402000. + 2400.*a['rounds'] + sum(a['idle_ns']),
                     benchmark_duration_ns=2000. + 2400.*a['rounds'] + sum(a['idle_ns'])) for a in group]
    copies = deepcopy(arms); result = fit_timing(None, copies, fake_compile)
    assert result['all_fit']
    for a in copies:
        assert max(a['idle_ns']) - min(a['idle_ns']) <= 4.
        assert abs(a['verified_duration_ns'] - TARGET_NS) <= TOLERANCE_NS
        assert a['verified_schedule_ns'] > 400_000.
    def too_long(_, group):
        return [dict(name=a['name'], total_schedule_ns=420_000., benchmark_duration_ns=20_000.) for a in group]
    copies = deepcopy(arms); assert not fit_timing(None, copies, too_long)['all_fit']
    assert all(not a['fits'] and sum(a['idle_ns']) == 0 for a in copies)


def check_mapping(arms, architecture):
    edges = {frozenset(locus) for impl in architecture['gates']['cz']['implementations'].values()
             for locus in impl['loci'] if len(locus) == 2}
    for a in arms:
        assert all(frozenset(i['locus']) in edges for i in instructions(a) if i['name'] == 'cz'), a['name']


def get_compiler(metadata):
    from dotenv import load_dotenv
    from iqm.pulla.pulla import Pulla
    from iqm.cpc.core.observation.observation_loading_rules import LatestFromStash
    load_dotenv(HERE.parent / '.env')
    if not os.environ.get('IQM_TOKEN'): raise RuntimeError('Set IQM_TOKEN in environment or .env')
    pulla = Pulla('https://resonance.iqm.tech', quantum_computer='emerald')
    stash = pulla.get_calibration_stash(UUID(metadata['calibration_set_id']))
    return pulla.get_standard_compiler(loading_rules=[LatestFromStash(stash)])


def analyze(output):
    plan = json.loads((output / 'plan.json').read_text()); summaries = {}
    for a in plan['arms']:
        files = sorted((output / a['name']).glob('block*.npz'))
        if len(files) != plan['metadata']['blocks']: raise RuntimeError(f'Incomplete data for {a["name"]}')
        failures = raw = shots = 0
        for path in files:
            with np.load(path) as saved:
                failures += int(np.any(saved['predictions'] != saved['obs_flips'], axis=1).sum())
                raw += int(saved['obs_flips'].any(axis=1).sum()); shots += len(saved['obs_flips'])
        probability = failures/shots; ci = interval(failures, shots)
        # Clamp floating-point roundoff at zero/all failures so plotting never
        # receives a tiny negative error-bar length.
        ci = [min(probability, ci[0]), max(probability, ci[1])]
        summaries[a['name']] = dict(shots=shots, failures=failures, p_error=probability,
                                   wilson95=ci, raw_probability=raw/shots)
    assert len({s['shots'] for s in summaries.values()}) == 1
    (output / 'results.json').write_text(json.dumps(dict(metadata=plan['metadata'], summary=summaries), indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, basis in zip(axes, ('Z', 'X')):
        rows = [summaries[f'{basis}_r{r}'] for r in range(4)]
        ps = np.array([s['p_error'] for s in rows]); cis = np.array([s['wilson95'] for s in rows])
        ax.errorbar(range(4), ps, yerr=np.maximum(0., np.stack([ps-cis[:, 0], cis[:, 1]-ps])), marker='o', capsize=4, label='Encoded; Wilson 95% CI')
        base = summaries[f'{basis}_single']; ax.axhline(base['p_error'], color='tab:orange', label='Unencoded qubit')
        ax.axhspan(*base['wilson95'], color='tab:orange', alpha=.15)
        suffix = 'incl. readout' if plan['timing']['includes_final_readout'] else 'before readout'
        ax.set(title=f'{basis}-memory, 15 μs after prep ({suffix})', xlabel='Intermediate extraction rounds', xticks=range(4))
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    axes[0].set_ylabel('Memory error probability'); fig.tight_layout()
    fig.savefig(output / 'ler_vs_rounds.png', dpi=180); plt.close(fig)
    print(json.dumps(summaries, indent=2))


def execute(output, plan, timeout):
    """Resumable native IQM pipeline; no correction gates, all history saved."""
    arms = plan['arms']; meta = plan['metadata']
    if plan['timing'].get('includes_initialization', True) or plan['timing'].get('includes_state_preparation', True):
        raise RuntimeError('Old total-duration plan: regenerate using the post-preparation timing definition')
    if plan['timing'].get('includes_final_readout', True):
        raise RuntimeError('This benchmark ends before final readout: regenerate the plan')
    if not plan['timing'].get('prefix_terminal_probes_removed', False):
        raise RuntimeError('Old prefix timing probe format: regenerate the plan')
    if not plan['timing']['all_fit']: raise RuntimeError('15 us infeasible: no hardware submission allowed')
    validate(arms)
    compiler = get_compiler(meta); rows = compile_durations(compiler, arms)
    (output / 'pre_submission_timing.json').write_text(json.dumps(rows, indent=2))
    if not all(-1e-6 <= TARGET_NS - row['benchmark_duration_ns'] <= TOLERANCE_NS for row in rows):
        raise RuntimeError('Recompiled schedule differs from target; no hardware submission allowed')
    from iqm.iqm_client import IQMClient
    client = IQMClient('https://resonance.iqm.tech', quantum_computer='emerald')
    dqa = client.get_dynamic_quantum_architecture(UUID(meta['calibration_set_id']))
    check_mapping(arms, dqa.model_dump(mode='json'))
    for block in range(meta['blocks']):
        ordered = arms[block % len(arms):] + arms[:block % len(arms)]
        raw_paths = [output / f'block{block}_{a["name"]}_raw.json' for a in ordered]
        if all(p.exists() for p in raw_paths): raws = [json.loads(p.read_text()) for p in raw_paths]
        else:
            job_path = output / f'block{block}_job.json'
            if job_path.exists():
                saved = json.loads(job_path.read_text()); assert saved['arms'] == [a['name'] for a in ordered]
                job = client.get_job(UUID(saved['job_id']))
            else:
                job = client.submit_circuits(native_circuits(ordered), shots=meta['shots_per_arm_per_block'], calibration_set_id=UUID(meta['calibration_set_id']))
                job_path.write_text(json.dumps(dict(job_id=str(job.job_id), arms=[a['name'] for a in ordered]), indent=2))
            raws = collect_job_results(job, len(ordered), timeout, output / f'block{block}_status.json')
            for path, raw in zip(raw_paths, raws): path.write_text(json.dumps(raw))
        for a, raw in zip(ordered, raws):
            directory = output / a['name']; ref = stim.Circuit.from_file(str(directory / 'reference.stim'))
            measurements, dets, obs = extract(raw, ref)
            assert len(obs) == meta['shots_per_arm_per_block']
            if a['kind'] == 'single': pred = np.zeros_like(obs)
            else:
                dem = stim.DetectorErrorModel.from_file(str(directory / 'model.dem'))
                matcher = pymatching.Matching.from_detector_error_model(dem, enable_correlations=meta['correlated'])
                pred = matcher.decode_batch(dets, enable_correlations=meta['correlated'])
            np.savez_compressed(directory / f'block{block}.npz', raw_meas=measurements, det_events=dets, obs_flips=obs, predictions=pred)
    analyze(output)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--self-test', action='store_true'); mode.add_argument('--prepare', action='store_true')
    mode.add_argument('--hardware', action='store_true')
    p.add_argument('--resume-dir', type=Path); p.add_argument('--source-run', type=Path, default=SOURCE)
    p.add_argument('--shots', type=int, default=1000); p.add_argument('--blocks', type=int, default=2)
    p.add_argument('--baseline-qubit'); p.add_argument('--idle-error-per-us', type=float, default=.003)
    p.add_argument('--correlated', action='store_true'); p.add_argument('--timeout-secs', type=float, default=600.)
    args = p.parse_args()
    if args.shots < 1 or args.blocks < 1 or not 0 <= args.idle_error_per_us < .75 or args.timeout_secs <= 0: p.error('Invalid shots/blocks/noise/timeout')
    if args.hardware:
        if not args.resume_dir: p.error('--hardware requires --resume-dir from --prepare')
        execute(args.resume_dir, json.loads((args.resume_dir / 'plan.json').read_text()), args.timeout_secs); return 0
    if args.resume_dir: p.error('--resume-dir is only for hardware')
    source = json.loads((args.source_run / 'plan.json').read_text()); previous = source['metadata']
    if previous['device'] != 'emerald': p.error('Source must be an Emerald experiment')
    layout = {int(k): v for k, v in previous['layout'].items()}
    arms = make_arms(layout, args.baseline_qubit)
    for a in arms: a['include_final_readout'] = False
    architecture = json.loads((args.source_run / 'architecture.json').read_text())
    check_mapping(arms, architecture); validate(arms); timing_self_test(arms)
    if args.self_test: return 0
    output = HERE.parent / 'output/memory_15us' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    output.mkdir(parents=True)
    meta = dict(device='emerald', calibration_set_id=previous['calibration_set_id'], layout=layout,
                shots_per_arm_per_block=args.shots, blocks=args.blocks, correlated=args.correlated,
                target_post_preparation_ns=TARGET_NS, includes_final_readout=False,
                idle_error_per_us=args.idle_error_per_us, noise=DEFAULT_NOISE,
                zero_swaps=True, no_midcircuit_reset=True, preparation='Common unitary encoded known-state preparation for every round count; non-fault-tolerant',
                decoder='All 8 checks at each round plus 4 final data-derived checks; no physical corrections',
                timing='From the global barrier AFTER encoding to the barrier BEFORE final basis rotation/readout; boundary times from compiled prefixes, no guessed reset duration',
                idle_model='Effective homogeneous Pauli depolarization per us on all active qubits; not calibrated T1/T2',
                server_timestamp_verified=False)
    timing = fit_timing(get_compiler(meta), arms)
    validate(arms)
    plan = dict(metadata=meta, timing=timing, arms=arms)
    (output / 'plan.json').write_text(json.dumps(plan, indent=2))
    report = ['# Fixed-duration quantum memory feasibility', '',
              'Target: **15000 ns AFTER initialization and state preparation, BEFORE final basis rotation/readout**.', '',
              'Sample-grid tolerance: at most 4 ns below target; no accepted overrun. Actual durations are listed explicitly.', '',
              '| Arm | Unpadded window (ns) | Verified window (ns) | Full schedule (ns) | Idle slots (ns) | Fits |', '|---|---:|---:|---:|---|---|']
    for a in arms:
        directory = output / a['name']; directory.mkdir()
        model(a).to_file(str(directory / 'reference.stim'))
        noisy = model(a, True, args.idle_error_per_us); noisy.to_file(str(directory / 'model.stim'))
        noisy.detector_error_model(decompose_errors=True).to_file(str(directory / 'model.dem'))
        (directory / 'native_instructions.json').write_text(json.dumps(instructions(a), indent=2))
        report.append(f'| {a["name"]} | {a["unpadded_duration_ns"]:.1f} | {a.get("verified_duration_ns", a["unpadded_duration_ns"]):.1f} | {a.get("verified_schedule_ns", a["unpadded_schedule_ns"]):.1f} | {a["idle_ns"]} | {a["fits"]} |')
    report += ['', 'Excluded boundary costs, from compiled prefixes:', '',
               '| Arm | Initialization (ns) | State preparation (ns) | Final rotation/readout (ns) |', '|---|---:|---:|---:|']
    for a in arms:
        if 'timing_details' in a:
            t = a['timing_details']
            report.append(f'| {a["name"]} | {t["initialization_ns"]:.1f} | {t["state_preparation_ns"]:.1f} | {t["final_readout_ns"]:.1f} |')
    if not timing['all_fit']:
        report += ['', '**Hardware submission blocked:** at least one post-preparation configuration cannot fit. No pulse implementation was modified. No LER exists for this unrun benchmark.']
    report += ['', 'Durations come from local PulLA schedules using the pinned calibration. They are not server execution timestamps.']
    (output / 'feasibility.md').write_text('\n'.join(report) + '\n')
    print('\n'.join(report), flush=True); print(f'Saved: {output}; no hardware submitted.', flush=True)
    return 0 if timing['all_fit'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
