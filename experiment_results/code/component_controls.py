"""Offline reconstruction and analysis of four full-schedule-matched controls.

op, instructions and make_arms are extracted from the validated original
component experiment, omitting only its fifth zero-round arm. No hardware calls.
"""

from copy import deepcopy
from collections import Counter
import json
import numpy as np
import pymatching
from scipy.stats import fisher_exact
import stim
from surface_code import make_stim_circuit
from utilities import extract

ARM_NAMES = ('r1_idle', 'r1_cz_echo', 'r1_ancilla_measure', 'r3')
LABELS = ('1 round + idle', '1 round + CZ echo', '1 round + ancilla measurement', '3 rounds')

def op(name,locus,**args):
    return dict(name=name,locus=list(locus),args=args)


def make_arms(source_plan):
    layout={int(k):v for k,v in source_plan['metadata']['layout'].items()}
    one=next(a for a in source_plan['arms'] if a['rounds']==1)
    three=next(a for a in source_plan['arms'] if a['rounds']==3)
    # Remove previous timing padding and collapse consecutive global barriers.
    bare=[]
    for item in one['instructions']:
        if item['name']=='delay': continue
        if item['name']=='barrier' and bare and bare[-1]['name']=='barrier': continue
        bare.append(deepcopy(item))
    boundary=max(i for i,item in enumerate(bare) if item['name']=='barrier')
    prefix,readout=bare[:boundary+1],bare[boundary+1:]
    assert len(readout)==9 and all(i['name']=='measure' for i in readout)
    all_q=list(bare[boundary]['locus'])
    assert set(all_q)==set(layout.values())
    source=make_stim_circuit(3,1,basis='Z',no_reset=True)
    coords=source.get_final_qubit_coordinates()
    ancillas={q for q,xy in coords.items() if not (int(xy[0])%2 and int(xy[1])%2)}
    layers=[]
    for inst in source.flattened():
        if inst.name=='CX':
            qs=[t.value for t in inst.targets_copy()]
            assert len(set(qs))==len(qs), 'Expected disjoint gates within each layer'
            layers.append([(qs[i],qs[i+1]) for i in range(0,len(qs),2)])
    assert sum(map(len,layers))==24
    echo_ops=[]; echo_pairs=[]; echo=stim.Circuit()
    for layer in layers:
        for repeat in range(2):
            for a,b in layer:
                echo_ops.append(op('cz',(layout[a],layout[b])))
                echo.append('CZ',[a,b]); echo_pairs.append([a,b])
            echo_ops.append(op('barrier',all_q))
    assert stim.Tableau.from_circuit(echo)==stim.Tableau(echo.num_qubits), 'CZ block must be exact identity'
    added_cz=sum(item['name']=='cz' for item in echo_ops)
    assert added_cz==sum(i['name']=='cz' for i in three['instructions'])-sum(i['name']=='cz' for i in bare)==48
    aux=[]
    for pass_index in range(2):
        for q in sorted(ancillas):
            aux.append(op('measure',(layout[q],),key=f'aux_{pass_index}_s{q}'))
        aux.append(op('barrier',all_q))
    first_measured={i['locus'][0] for i in prefix if i['name']=='measure'}
    assert first_measured=={layout[q] for q in ancillas}
    assert all(i['name']=='barrier' or (i['name']=='measure' and set(i['locus'])<=first_measured) for i in aux)
    arms=[dict(name='r1_idle',rounds=1,component='none',base_ops=prefix+readout),
          dict(name='r1_cz_echo',rounds=1,component='48_CZ_identity_echo',base_ops=prefix+echo_ops+readout,echo_pairs=echo_pairs),
          dict(name='r1_ancilla_measure',rounds=1,component='16_extra_ancilla_Z_readouts',base_ops=prefix+aux+readout),
          dict(name='r3',rounds=3,component='three_rounds',base_ops=deepcopy(three['instructions']))]
    for arm in arms:
        i=max(j for j,item in enumerate(arm['base_ops']) if item['name']=='barrier')
        assert len(arm['base_ops'][i+1:])==9
        arm['padding_ns']=0.
        arm['boundary']=i
    return arms,layout


def instructions(arm):
    items=deepcopy(arm['base_ops']); i=arm['boundary']; padding=arm['padding_ns']
    if padding>0:
        items[i+1:i+1]=[op('delay',items[i]['locus'],duration=float(padding*1e-9)),op('barrier',items[i]['locus'])]
    return items


def analyze_controls(root, output, edges, calibration_id, wilson):
    run = root/'data/component_controls'
    plan = json.loads((run/'plan.json').read_text()); meta = plan['metadata']
    saved = json.loads((run/'results.json').read_text())['summary']
    source = json.loads((root/'data/source_patch/plan.json').read_text())
    rebuilt, layout = make_arms(source)
    expected = {a['name']: a for a in rebuilt}
    arms = {a['name']: a for a in plan['arms'] if a['name'] in ARM_NAMES}
    assert len(arms) == 4 and set(expected) == set(arms)
    assert meta['layout'] == source['metadata']['layout'] and meta['calibration_set_id'] == calibration_id
    assert meta['timing']['matched_locally'] and not meta['timing']['server_verified']
    timing = {t['arm']: t for t in meta['timing']['final']}
    record_mode = meta['timing'].get('matching_basis') == 'first_round_record_to_final_readout_boundary'
    target = meta['timing'].get('target_record_window_ns') if record_mode else meta['timing']['target_full_schedule_ns']
    job_ids_by_arm = {}
    for name in ARM_NAMES:
        acquisition = run/plan.get('acquisition_by_arm', {}).get(name, '.')
        job_ids_by_arm[name] = []
        for block in range(meta['blocks']):
            job = json.loads((acquisition/f'block{block}_job.json').read_text())
            status = json.loads((acquisition/f'block{block}_status.json').read_text())
            assert status['job_id'] == job['job_id'] and status['status'] == 'completed'
            assert name in job['arms']
            job_ids_by_arm[name].append(job['job_id'])
    summaries = []; timings = []; details = {}; native_by_arm = {}; models = {}
    for name in ARM_NAMES:
        a = arms[name]; folder = run/name
        assert a['base_ops'] == expected[name]['base_ops'] and a['boundary'] == expected[name]['boundary']
        ops = instructions(a); native_by_arm[name] = ops
        assert all(frozenset(i['locus']) in edges for i in ops if i['name'] == 'cz')
        assert all(i['name'] in ('prx', 'cz', 'barrier', 'delay', 'measure') for i in ops)
        counts = Counter(i['name'] for i in ops)
        assert counts['cz'] == timing[name]['compiled_cz_count'] == (67 if name in ('r1_cz_echo', 'r3') else 19)
        checked_duration = timing[name]['record_window_ns'] if record_mode else timing[name]['full_schedule_ns']
        assert abs(checked_duration-target) <= meta['timing']['tolerance_ns']
        reference = stim.Circuit.from_file(str(folder/'reference.stim'))
        model = stim.Circuit.from_file(str(folder/'model.stim'))
        dem = stim.DetectorErrorModel.from_file(str(folder/'model.dem')); models[name] = dem
        assert reference.without_noise() == model.without_noise()
        assert reference.num_measurements == 8*a['rounds']+9 and reference.num_detectors == 8*a['rounds']
        main_keys = [i['args']['key'] for i in ops if i['name'] == 'measure' and i['args']['key'].startswith('m_')]
        aux_keys = [i['args']['key'] for i in ops if i['name'] == 'measure' and i['args']['key'].startswith('aux_')]
        # The earlier optimizer reorders commuting measurements. Named keys,
        # rather than physical instruction order, define Stim record columns.
        assert len(main_keys) == reference.num_measurements
        assert set(main_keys) == {f'm_{j:03d}' for j in range(reference.num_measurements)}
        assert len(aux_keys) == (16 if name == 'r1_ancilla_measure' else 0)
        matcher = pymatching.Matching.from_detector_error_model(dem, enable_correlations=meta['correlated'])
        n = failures = raw_flips = 0; blocks = []
        for block in range(meta['blocks']):
            with np.load(folder/f'block{block}.npz', allow_pickle=False) as data:
                ds, os_ = reference.compile_m2d_converter().convert(measurements=data['raw_meas'].astype(bool), separate_observables=True)
                assert np.array_equal(ds, data['det_events']) and np.array_equal(os_, data['obs_flips'])
                prediction = matcher.decode_batch(ds, enable_correlations=meta['correlated'])
                assert np.array_equal(prediction, data['predictions'])
                assert len(os_) == meta['shots_per_arm_per_block']
                if aux_keys:
                    raw = json.loads((run/f'block{block}_r1_ancilla_measure_raw.json').read_text())
                    assert set(raw) == set(main_keys+aux_keys)
                    assert all(np.asarray(raw[k]).shape == (len(os_), 1) for k in aux_keys)
                    measured, rd, ro = extract({k:raw[k] for k in main_keys}, reference)
                    assert np.array_equal(measured, data['raw_meas']) and np.array_equal(rd, ds) and np.array_equal(ro, os_)
                bad = int(np.any(prediction != os_, axis=1).sum())
                raw_bad = int(os_.any(axis=1).sum())
                blocks.append(dict(block=block, shots=len(os_), failures=bad, raw_flips=raw_bad))
                n += len(os_); failures += bad; raw_flips += raw_bad
        s = saved[name]
        assert (n, failures, raw_flips) == (s['shots'], s['failures'], s['raw_flips'])
        assert failures/n == s['p_error'] and raw_flips/n == s['raw_probability']
        assert np.allclose(wilson(failures, n), s['wilson95'], atol=1e-12, rtol=0)
        summaries.append(dict(dataset='component_controls', arm=name, basis='Z', kind='legacy_Z_memory_component_control', rounds=a['rounds'],
            shots=n, failures=failures, ler=s['p_error'], wilson95_low=s['wilson95'][0], wilson95_high=s['wilson95'][1],
            raw_flips=raw_flips, raw_probability=s['raw_probability'], calibration_set_id=meta['calibration_set_id'], job_ids='|'.join(job_ids_by_arm[name])))
        timings.append(dict(dataset='component_controls', arm=name, requested_memory_ns=target if record_mode else '', compiled_memory_ns=checked_duration if record_mode else '',
            initialization_ns='', state_preparation_ns='', final_readout_ns='', full_schedule_ns=timing[name]['full_schedule_ns'],
            idle_slots_ns=json.dumps([a['padding_ns']]), native_prx=counts['prx'], native_cz=counts['cz'], native_measurements=counts['measure'],
            target_full_schedule_ns='' if record_mode else target, timing_definition='First complete syndrome record to pre-readout barrier; NOT a shared independent encoded-state boundary' if record_mode else 'Full schedule matched; separate memory boundaries not archived'))
        details[name] = dict(blocks=blocks, counts=dict(counts), padding_ns=a['padding_ns'], final_schedule_ns=timing[name]['full_schedule_ns'])
    assert models['r1_idle'].approx_equals(models['r1_cz_echo'], atol=1e-12)
    assert models['r1_idle'].approx_equals(models['r1_ancilla_measure'], atol=1e-12)
    assert details['r1_cz_echo']['counts']['cz']-details['r1_idle']['counts']['cz'] == 48
    assert details['r1_ancilla_measure']['counts']['measure']-details['r1_idle']['counts']['measure'] == 16

    base = summaries[0]; contrasts = []
    for s in summaries[1:]:
        delta = s['ler']-base['ler']
        lo = delta-np.hypot(s['ler']-s['wilson95_low'], base['wilson95_high']-base['ler'])
        hi = delta+np.hypot(s['wilson95_high']-s['ler'], base['ler']-base['wilson95_low'])
        p = float(fisher_exact([[s['failures'], s['shots']-s['failures']], [base['failures'], base['shots']-base['failures']]]).pvalue)
        contrasts.append(dict(arm=s['arm'], delta_vs_idle=delta, difference95_low=float(lo), difference95_high=float(hi), fisher_p=p))
    running = 0.
    for i, c in enumerate(sorted(contrasts, key=lambda c:c['fisher_p'])):
        running = max(running, min(1., (len(contrasts)-i)*c['fisher_p'])); c['holm_p_3_contrasts'] = running
    report = dict(arms=list(ARM_NAMES), shots=sum(s['shots'] for s in summaries),
        details=details, contrasts=contrasts, original_primary_decoder_preserved=True, three_one_round_dems_identical=True,
        auxiliary_measurements_preserved_but_excluded_from_decoder=True, fixed_patch_recipe_reproduced=True,
        timing_scope=f'First complete syndrome record to final readout barrier, target {target/1000:g} us; first three arms rerun, r3 is historical; NOT a common independent encoding boundary' if record_mode else 'Full local schedules ~410.872 us, NOT the 15 us memory benchmark; no isolated memory-window duration was archived',
        contrast_scope='Derived three displayed contrasts only, unpaired shots, Fisher test + Holm over 3; original archive retains its original statistics unchanged')
    (output/'tables/component_controls.json').write_text(json.dumps(dict(summary=summaries, timing=timings, analysis=report), indent=2)+'\n')
    (output/'tables/component_native_instructions.json').write_text(json.dumps(native_by_arm, indent=2)+'\n')
    lines = ['# Four component controls: '+('first-syndrome-record window matching' if record_mode else 'local full-schedule matching'), '',
        'Z-memory, distinct from the common unitary-encoder 15 us benchmark. Measured counts, primary DEMs and acquisition metadata are preserved.',
        report['timing_scope'], '',
        '| Arm | Failures / shots | LER (Wilson 95% CI) | Raw parity error | CZ | Padding | Full schedule |',
        '|---|---:|---|---:|---:|---:|---:|']
    for label, s in zip(LABELS, summaries):
        info = details[s['arm']]
        lines.append(f'| {label} | {s["failures"]}/{s["shots"]} | {s["ler"]:.2%} [{s["wilson95_low"]:.2%}, {s["wilson95_high"]:.2%}] | {s["raw_probability"]:.2%} | {info["counts"]["cz"]} | {info["padding_ns"]/1000:g} us | {info["final_schedule_ns"]/1000:g} us |')
    lines += ['', '| Contrast vs 1 round + idle | LER difference (pp) | Unadjusted 95% difference CI (pp) | Holm p (3 contrasts) |', '|---|---:|---|---:|']
    for c in contrasts:
        lines.append(f'| {c["arm"]} | {100*c["delta_vs_idle"]:+.2f} | [{100*c["difference95_low"]:+.2f}, {100*c["difference95_high"]:+.2f}] | {c["holm_p_3_contrasts"]:.4g} |')
    lines += ['', '## Interpretation', '',
        f'- Three-round minus one-round-idle LER difference: {100*next(c["delta_vs_idle"] for c in contrasts if c["arm"] == "r3"):+.2f} percentage points. Refer to the intervals and adjusted tests above.',
        '- A non-significant contrast is not an equivalence test or evidence that gate/readout noise is absent.',
        '- The CZ intervention is 48 extra gates in ideal CZ-CZ identity pairs on the 24 check links, not arbitrary unpaired entangling gates.',
        '- The measurement intervention repeats the already-measured ancillas twice in Z, without re-entangling data or adding H. It records 16 auxiliary bits which are deliberately excluded from the primary decoder.',
        '- All three one-round controls have the identical primary DEM; three rounds use their own full-history DEM.',
        '- Matching uses the first complete syndrome record, not a common independently encoded-state barrier. The reused three-round measurements come from earlier jobs; drift is a confound.' if record_mode else '- Full schedule durations, including initialization/readout, differ by at most 4 ns. Operation additions replace waiting; the contrast measures their net effect, not isolated error per gate.',
        '- These interventions do not reproduce the full gate/measurement/basis-rotation sequence of repeated extraction, so the gap cannot be uniquely attributed to CZ or measurement alone.',
        '- Only the latest three one-round measurements and the required historical three-round reference are retained here. Superseded trials remain in diagnostics.',
        '- A record-availability timestamp is not a separately established logical-state preparation boundary. Server timestamps and full waveforms remain unavailable.' if record_mode else '- Missing isolated memory-boundary timings are blank in timing.csv; they are not inferred by subtracting an assumed reset time. Server timestamps and full waveforms remain unavailable.']
    (output/'COMPONENT_CONTROLS.md').write_text('\n'.join(lines)+'\n')
    plot_controls(summaries, contrasts, output/'plots/component_controls.png', record_mode)
    return summaries, timings, report


def plot_controls(rows, contrasts, path, record_mode=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    p = np.array([s['ler'] for s in rows])*100
    low = np.array([s['wilson95_low'] for s in rows])*100
    high = np.array([s['wilson95_high'] for s in rows])*100
    x = np.arange(4)
    axes[0].errorbar(x, p, yerr=np.stack([p-low, high-p]), fmt='o', capsize=5, label='Decoded, Wilson 95% CI')
    axes[0].plot(x, [s['raw_probability']*100 for s in rows], 's--', label='Raw logical parity error', alpha=.8)
    axes[0].set(xticks=x, xticklabels=['1 round\n+ idle', '1 round\n+ CZ echo', '1 round\n+ ancilla reads', '3 rounds'],
                ylabel='Error probability (%)', title='Z-memory, matched first-record windows' if record_mode else 'Z-memory, same local full-schedule duration', ylim=(0, 60))
    axes[0].legend(fontsize=8); axes[0].grid(alpha=.2)
    d = np.array([c['delta_vs_idle'] for c in contrasts])*100
    lo = np.array([c['difference95_low'] for c in contrasts])*100
    hi = np.array([c['difference95_high'] for c in contrasts])*100
    axes[1].errorbar(range(3), d, yerr=np.stack([d-lo, hi-d]), fmt='o', capsize=5)
    axes[1].axhline(0, color='gray', linestyle='--')
    axes[1].set(xticks=range(3), xticklabels=['+ CZ echo', '+ ancilla reads', '3 rounds'],
                ylabel='LER difference vs 1 round + idle (pp)', title='Contrasts with 95% difference intervals')
    axes[1].grid(alpha=.2)
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)
