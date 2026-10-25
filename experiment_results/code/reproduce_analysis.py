"""Verify and analyze the archived measurements, strictly offline.

Run from any directory: python /path/to/submission/code/reproduce_analysis.py
Writes derived tables/plots/reports only; never edits data/ or calls IQM.
No new random samples, simulations, hardware jobs or parameter fitting.
"""

import argparse
import ast
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pymatching
import stim

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ('memory_15us', 'memory_near_zero', 'memory_5us')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def wilson(k, n):
    p = k/n; z = 1.959963984540054; den = 1+z*z/n
    center = (p+z*z/(2*n))/den
    half = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [max(0., min(p, center-half)), min(1., max(p, center+half))]


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT)
    args = parser.parse_args(); output = args.output_dir.resolve()
    # Outputs must never overwrite the immutable archival measurement tree.
    if output == (ROOT/'data') or (ROOT/'data') in output.parents:
        parser.error('Choose an output directory outside data/')
    for path in (output/'tables', output/'plots'): path.mkdir(parents=True, exist_ok=True)
    provenance = json.loads((ROOT/'PROVENANCE.json').read_text())
    copies = [r for r in provenance['files'] if r['mode'] == 'byte-identical copy']
    for item in copies:
        assert sha(ROOT/item['destination']) == item['sha256'], item['destination']
    for path in (ROOT/'code').glob('*.py'): ast.parse(path.read_text())

    # Import only submission-local definitions. This exact Clifford construction
    # checks encoded stabilizers and bus cleanup; it generates no shot samples.
    import memory_benchmark as benchmark
    source = json.loads((ROOT/'data/source_patch/plan.json').read_text())
    layout = {int(q): v for q, v in source['metadata']['layout'].items()}
    assert len(layout) == len(set(layout.values())) == 17
    architecture = json.loads((ROOT/'data/source_patch/architecture.json').read_text())
    assert architecture['calibration_set_id'] == source['metadata']['calibration_set_id']
    edges = {frozenset(locus) for impl in architecture['gates']['cz']['implementations'].values()
             for locus in impl['loci'] if len(locus) == 2}
    baseline_plan = json.loads((ROOT/'data/memory_15us/plan.json').read_text())
    single = next(a for a in baseline_plan['arms'] if a['name'] == 'Z_single')
    recipes = {a['name']: a for a in benchmark.make_arms(layout, layout[single['data_qubit']])}

    summaries = []; timings = []; native_by_arm = {}; named_timings = {}; total_shots = 0
    for dataset in DATASETS:
        run = ROOT/'data'/dataset
        plan = json.loads((run/'plan.json').read_text()); meta = plan['metadata']
        saved = json.loads((run/'results.json').read_text())['summary']
        rows = {r['name']: r for r in json.loads((run/'pre_submission_timing.json').read_text())}
        assert meta['layout'] == baseline_plan['metadata']['layout']
        assert meta['calibration_set_id'] == architecture['calibration_set_id']
        assert meta['device'] == 'emerald' and not meta['includes_final_readout']
        job_ids = []
        for b in range(meta['blocks']):
            job = json.loads((run/f'block{b}_job.json').read_text())
            status = json.loads((run/f'block{b}_status.json').read_text())
            assert job['job_id'] == status['job_id'] and status['status'] == 'completed'
            assert set(job['arms']) == {a['name'] for a in plan['arms']}
            job_ids.append(job['job_id'])
        for arm in plan['arms']:
            name = arm['name']; folder = run/name
            reference = stim.Circuit.from_file(str(folder/'reference.stim'))
            dem = stim.DetectorErrorModel.from_file(str(folder/'model.dem'))
            model = stim.Circuit.from_file(str(folder/'model.stim'))
            assert reference.without_noise() == model.without_noise()
            native = json.loads((folder/'native_instructions.json').read_text())
            assert benchmark.instructions(arm) == native
            recipe = recipes[arm.get('source_arm', name)]
            for key in ('prep_ops', 'round_ops', 'readout_ops', 'active_locus', 'ancilla_locus'):
                assert arm[key] == recipe[key], (dataset, name, key)
            reconstructed = dict(recipe, idle_ns=arm['idle_ns'])
            assert benchmark.instructions(reconstructed) == native
            assert all(frozenset(i['locus']) in edges for i in native if i['name'] == 'cz')
            assert not any(i['name'] in ('swap', 'move', 'reset', 'r') for i in native)
            counts = Counter(i['name'] for i in native)
            expected_cz = 16+24*arm['rounds'] if arm['kind'] == 'encoded' else 0
            assert counts['cz'] == expected_cz
            assert reference.num_measurements == (9+8*arm['rounds'] if arm['kind'] == 'encoded' else 1)
            assert reference.num_detectors == (4+8*arm['rounds'] if arm['kind'] == 'encoded' else 0)
            measured_keys = [i['args']['key'] for i in native if i['name'] == 'measure']
            assert measured_keys == [f'm_{j:03d}' for j in range(reference.num_measurements)]
            matcher = None if arm['kind'] == 'single' else pymatching.Matching.from_detector_error_model(dem, enable_correlations=meta['correlated'])
            n = failures = raw_flips = 0
            for block in range(meta['blocks']):
                with np.load(folder/f'block{block}.npz', allow_pickle=False) as data:
                    measurements = data['raw_meas'].astype(bool)
                    ds, os_ = reference.compile_m2d_converter().convert(measurements=measurements, separate_observables=True)
                    assert np.array_equal(ds, data['det_events']) and np.array_equal(os_, data['obs_flips'])
                    prediction = np.zeros_like(os_) if matcher is None else matcher.decode_batch(ds, enable_correlations=meta['correlated'])
                    assert np.array_equal(prediction, data['predictions']), (dataset, name, block)
                    assert len(os_) == meta['shots_per_arm_per_block']
                    n += len(os_); failures += int(np.any(prediction != os_, axis=1).sum())
                    raw_flips += int(os_.any(axis=1).sum())
            s = saved[name]
            assert (n, failures) == (s['shots'], s['failures'])
            assert failures/n == s['p_error'] and raw_flips/n == s['raw_probability']
            assert np.allclose(wilson(failures, n), s['wilson95'], atol=1e-12, rtol=0)
            total_shots += n
            t = rows[name]
            assert np.isclose(t['total_schedule_ns'], sum(t[k] for k in ('initialization_ns', 'state_preparation_ns', 'storage_before_readout_ns', 'final_readout_ns')), atol=1e-6, rtol=0)
            assert t['final_readout_ns'] > 0 and not t['includes_final_readout']
            requested = {'memory_15us': 15000., 'memory_near_zero': 0., 'memory_5us': 5000.}[dataset]
            actual = t['storage_before_readout_ns']
            assert (-1e-6 <= requested-actual <= 4.000001) if requested else (-1e-6 <= actual <= 8.)
            assert max(arm['idle_ns'])-min(arm['idle_ns']) <= 4.
            summaries.append(dict(dataset=dataset, arm=name, basis=arm['basis'], kind=arm['kind'], rounds=arm['rounds'],
                shots=n, failures=failures, ler=s['p_error'], wilson95_low=s['wilson95'][0], wilson95_high=s['wilson95'][1],
                raw_flips=raw_flips, raw_probability=s['raw_probability'], calibration_set_id=meta['calibration_set_id'],
                job_ids='|'.join(job_ids)))
            timings.append(dict(dataset=dataset, arm=name, requested_memory_ns=requested, compiled_memory_ns=actual,
                initialization_ns=t['initialization_ns'], state_preparation_ns=t['state_preparation_ns'], final_readout_ns=t['final_readout_ns'],
                full_schedule_ns=t['total_schedule_ns'], idle_slots_ns=json.dumps(arm['idle_ns']), native_prx=counts['prx'], native_cz=counts['cz'], native_measurements=counts['measure'],
                target_full_schedule_ns='', timing_definition='Post-preparation memory window, before final basis rotation/readout'))
            native_by_arm[(dataset, name)] = native; named_timings[(dataset, name)] = t

    x5 = native_by_arm[('memory_5us', 'X_r0_5us')]; x15 = native_by_arm[('memory_15us', 'X_r0')]
    diff = [j for j, (a, b) in enumerate(zip(x5, x15)) if a != b]
    assert len(x5) == len(x15) == 136 and diff == [107]
    assert x5[107]['name'] == x15[107]['name'] == 'delay'
    assert x5[107]['locus'] == x15[107]['locus']
    assert x5[107]['args']['duration'] == 5e-6 and x15[107]['args']['duration'] == 15e-6
    tx5 = named_timings[('memory_5us', 'X_r0_5us')]; tx15 = named_timings[('memory_15us', 'X_r0')]
    for key in ('initialization_ns', 'state_preparation_ns', 'final_readout_ns'):
        assert tx5[key] == tx15[key]
    assert tx15['total_schedule_ns']-tx5['total_schedule_ns'] == 10000.
    assert total_shots == 28000 and len(summaries) == 14

    from component_controls import analyze_controls
    control_rows, control_timings, control_report = analyze_controls(ROOT, output, edges, architecture['calibration_set_id'], wilson)
    summaries.extend(control_rows); timings.extend(control_timings); total_shots += control_report['shots']
    assert total_shots == 36000 and len(summaries) == 18

    write_csv(output/'tables/summary.csv', summaries); write_csv(output/'tables/timing.csv', timings)
    (output/'tables/summary.json').write_text(json.dumps(dict(summary=summaries, timing=timings), indent=2)+'\n')
    validation = dict(archived_copy_hashes_verified=len(copies), configurations=len(summaries), measured_shots_redecoded=total_shots,
        raw_to_detector_observable_alignment=True, original_mwpm_predictions_reproduced=True, original_failure_counts_and_intervals_reproduced=True,
        all_native_cz_links_valid=True, primary_14_native_sequences_equal_to_direct_no_routing_recipe=True,
        component_control_fixed_patch_recipe_reproduced=True, component_control_primary_decoders_reproduced=True,
        component_control_auxiliary_readouts_preserved=True,
        saved_timing_decompositions_consistent=True, x_5us_15us_native_difference_is_only_idle=True,
        complete_pulse_waveforms_available=False, server_schedule_equivalence_verified=False,
        network_calls=0, hardware_jobs_submitted=0, new_experiments=0)
    (output/'validation.json').write_text(json.dumps(validation, indent=2)+'\n')
    plot(summaries, output/'plots')
    lines = ['# Reproduced measured results', '',
        'All counts below are recomputed from the archived shots using their original DEM and MWPM settings. No tuning or new experiments.', '',
        '| Dataset | Arm | Failures / shots | LER | Wilson 95% CI |', '|---|---|---:|---:|---|']
    for s in summaries:
        lines.append(f'| {s["dataset"]} | {s["arm"]} | {s["failures"]}/{s["shots"]} | {s["ler"]:.2%} | [{s["wilson95_low"]:.2%}, {s["wilson95_high"]:.2%}] |')
    lines += ['', 'Timing values in tables/timing.csv are archived local pulse-compiler summaries, not server timestamps or full waveforms.',
              'The primary 14 native gate sequences match the direct no-routing recipe; the four earlier controls reproduce their saved fixed-patch intervention recipe. All CZ links exist in the saved architecture. This is gate-level verification.',
              'Component-control timing scope: '+control_report['timing_scope']+'. See COMPONENT_CONTROLS.md for details.',
              'Near-zero/5 us/15 us results come from separate jobs. Calibration ID is shared; drift is not ruled out.']
    (output/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(validation, indent=2)); print(f'Output: {output}')


def plot(rows, folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def errors(selected):
        p = np.array([s['ler'] for s in selected])
        low = np.array([s['wilson95_low'] for s in selected]); high = np.array([s['wilson95_high'] for s in selected])
        return p, np.maximum(0., np.stack([p-low, high-p]))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, basis in zip(axes, ('Z', 'X')):
        group = sorted((s for s in rows if s['dataset'] == 'memory_15us' and s['basis'] == basis and s['kind'] == 'encoded'), key=lambda s:s['rounds'])
        ps, error = errors(group)
        ax.errorbar(range(4), ps, yerr=error, marker='o', capsize=4, label='Encoded, Wilson 95% CI')
        base = next(s for s in rows if s['dataset'] == 'memory_15us' and s['basis'] == basis and s['kind'] == 'single')
        ax.axhline(base['ler'], color='tab:orange', label='Single-qubit baseline')
        ax.axhspan(base['wilson95_low'], base['wilson95_high'], color='tab:orange', alpha=.15)
        ax.set(title=f'{basis}-memory, target 15 μs', xlabel='Intermediate syndrome rounds', xticks=range(4), ylim=(0, .6))
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    axes[0].set_ylabel('Memory error probability')
    fig.tight_layout(); fig.savefig(folder/'memory_15us.png', dpi=180); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, basis in zip(axes, ('Z', 'X')):
        selected = [next(s for s in rows if s['dataset'] == dataset and s['basis'] == basis and s['kind'] == 'encoded' and s['rounds'] == 0)
                    for dataset in ('memory_near_zero', 'memory_5us', 'memory_15us')]
        ps, error = errors(selected)
        ax.errorbar([0, 5, 15], ps, yerr=error, marker='o', capsize=4)
        ax.set(title=f'{basis}-memory, zero extraction rounds', xlabel='Compiled storage window (μs)', xticks=[0, 5, 15], ylim=(0, .6))
        ax.grid(alpha=.2)
    axes[0].set_ylabel('Logical error probability')
    fig.suptitle('Identical preparation/readout; separate hardware jobs', fontsize=11)
    fig.tight_layout(); fig.savefig(folder/'zero_round_storage.png', dpi=180); plt.close(fig)


if __name__ == '__main__': main()
