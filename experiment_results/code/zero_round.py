"""Zero-round memory with no explicit storage delay, retaining the full encoder.

python code/zero_round.py --self-test
python code/zero_round.py --prepare
python code/zero_round.py --hardware --resume-dir DIR
python code/reproduce_analysis.py

Only --hardware submits. Preparation reads calibration and pulse-compiles.
This measures a preparation/readout/decoder error floor, not preparation alone.
Existing benchmark and pulse implementations are unchanged.
"""

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from uuid import UUID

import numpy as np
import pymatching
from scipy.stats import fisher_exact
import stim
import memory_benchmark as memory

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'data/memory_15us'
ARM_NAMES = ('Z_r0', 'X_r0')
CORE_KEYS = ('basis', 'kind', 'rounds', 'data_qubit', 'active_locus', 'prep_ops',
             'round_ops', 'ancilla_locus', 'readout_ops')


def signature(ops):
    return hashlib.sha256(json.dumps(ops, sort_keys=True).encode()).hexdigest()


def make_cases(source, storage_ns=0.):
    originals = {a['name']: a for a in source['arms']}
    cases = []
    for name in ARM_NAMES:
        original = originals[name]
        assert original['rounds'] == 0
        a = {key: deepcopy(original[key]) for key in CORE_KEYS}
        suffix = 'near0' if storage_ns == 0 else f'{storage_ns/1000:g}us'
        a.update(name=f'{name}_{suffix}', source_arm=name, idle_ns=[float(storage_ns)],
                 target_storage_ns=float(storage_ns),
                 include_final_readout=False,
                 preparation_sha256=signature(original['prep_ops']),
                 readout_sha256=signature(original['readout_ops']))
        assert signature(a['prep_ops']) == signature(original['prep_ops'])
        assert signature(a['readout_ops']) == signature(original['readout_ops'])
        # Only the storage delay changes; preserve every other operation.
        expected = [i for i in memory.instructions(original) if i['name'] != 'delay']
        assert [i for i in memory.instructions(a) if i['name'] != 'delay'] == expected
        assert sum(i['name'] == 'measure' for i in memory.instructions(a)) == 9
        cases.append(a)
    return cases


def validate(cases, layout):
    if len(cases) != 2 or {a['source_arm'] for a in cases} != set(ARM_NAMES):
        raise RuntimeError('This diagnostic requires exactly the two encoded Z/X zero-round arms')
    for a in cases:
        target = a.get('target_storage_ns', 0.)
        assert target >= 0 and a['kind'] == 'encoded' and a['rounds'] == 0 and a['idle_ns'] == [target]
        memory.native_preparation_check(a, layout)
    assert len({a.get('target_storage_ns', 0.) for a in cases}) == 1
    memory.validate(cases)


def verify_timing(cases, rows, max_near_zero_ns):
    assert len(cases) == len(rows)
    for a, row in zip(cases, rows):
        assert row['name'] == a['name']
        t = row['storage_before_readout_ns']
        target = a.get('target_storage_ns', 0.)
        if target == 0:
            if not -1e-6 <= t <= max_near_zero_ns:
                raise RuntimeError(f'{a["name"]}: near-zero interval is {t} ns, exceeds {max_near_zero_ns} ns')
        elif not -1e-6 <= target - t <= memory.TOLERANCE_NS:
            raise RuntimeError(f'{a["name"]}: compiled storage {t} ns differs from target {target} ns')
        assert row['final_readout_ns'] > 0


def score_interval(k, n):
    p = k/n; lo, hi = memory.interval(k, n)
    return [min(p, lo), max(p, hi)]


def analyze(output):
    plan = json.loads((output/'plan.json').read_text()); meta = plan['metadata']
    summaries = {}
    for a in plan['arms']:
        files = sorted((output/a['name']).glob('block*.npz'))
        if len(files) != meta['blocks']: raise RuntimeError(f'Incomplete blocks for {a["name"]}')
        total = dict(shots=0, failures=0, raw_flips=0, rescued_by_decoder=0, introduced_by_decoder=0)
        detector_counts = None; per_block = []
        for file in files:
            with np.load(file) as data:
                raw_bad = data['obs_flips'].any(axis=1)
                bad = np.any(data['predictions'] != data['obs_flips'], axis=1)
                n = len(bad)
                if n != meta['shots_per_arm_per_block']: raise RuntimeError('Shot count changed')
                total['shots'] += n; total['failures'] += int(bad.sum()); total['raw_flips'] += int(raw_bad.sum())
                total['rescued_by_decoder'] += int((raw_bad & ~bad).sum())
                total['introduced_by_decoder'] += int((~raw_bad & bad).sum())
                count = data['det_events'].sum(axis=0)
                detector_counts = count if detector_counts is None else detector_counts + count
                per_block.append(dict(block=file.stem, shots=n, failures=int(bad.sum()), raw_flips=int(raw_bad.sum())))
        total.update(p_error=total['failures']/total['shots'], wilson95=score_interval(total['failures'], total['shots']),
                     raw_probability=total['raw_flips']/total['shots'], raw_wilson95=score_interval(total['raw_flips'], total['shots']),
                     detector_rates=(detector_counts/total['shots']).tolist(), per_block=per_block)
        summaries[a['name']] = total
    assert len({s['shots'] for s in summaries.values()}) == 1
    comparisons = []
    for base in ARM_NAMES:
        arm = next(a for a in plan['arms'] if a['source_arm'] == base)
        near = summaries[arm['name']]
        control = summaries.get(f'{base}_15us')
        provenance = 'same diagnostic experiment, separate shots'
        if control is None:
            control = meta.get('historical_15us_summary', {}).get(base)
            provenance = 'historical job: drift and run-to-run differences may confound the comparison'
        if control:
            delta = near['p_error'] - control['p_error']; lo1, hi1 = near['wilson95']; lo2, hi2 = control['wilson95']
            ci = [delta - np.hypot(near['p_error']-lo1, hi2-control['p_error']),
                  delta + np.hypot(hi1-near['p_error'], control['p_error']-lo2)]
            pv = float(fisher_exact([[near['failures'], near['shots']-near['failures']],
                                    [control['failures'], control['shots']-control['failures']]]).pvalue)
            comparisons.append(dict(source_arm=base, diagnostic_minus_15us=delta, difference_ci95=ci,
                                    fisher_p_unadjusted=pv, control=control, provenance=provenance))
            if meta.get('target_post_preparation_ns', 0.) == 0:
                comparisons[-1]['near0_minus_15us'] = delta  # Preserve the existing zero-delay result field.
    result = dict(metadata=meta, timing=plan['timing'], summary=summaries, comparisons=comparisons)
    (output/'results.json').write_text(json.dumps(result, indent=2))
    storage_us = meta.get('target_post_preparation_ns', 0.) / 1000
    lines = [f'# Zero-round memory diagnostic: {storage_us:g} us storage', '',
             '| Arm | Shots | Decoded error | Wilson 95% CI | Raw parity error | Decoder rescued / introduced |',
             '|---|---:|---:|---|---:|---:|']
    for name, s in summaries.items():
        lo, hi = s['wilson95']
        lines.append(f'| {name} | {s["shots"]} | {s["p_error"]:.2%} | [{lo:.2%}, {hi:.2%}] | {s["raw_probability"]:.2%} | {s["rescued_by_decoder"]} / {s["introduced_by_decoder"]} |')
    lines += ['', '## Interpretation', '']
    if storage_us == 0:
        lines += ['- Low near-zero error bounds the combined preparation/readout/decoder error floor for the measured logical basis. It does not isolate preparation fidelity.',
                  '- High near-zero error means the long storage delay is not necessary for failure. Investigate encoding, initialization, readout and decoding before attributing the entire error to storage.']
    else:
        lines += [f'- This includes {storage_us:g} us of storage as well as preparation/readout errors. Compare with the separate near-zero experiment to assess the added storage damage.',
                  '- Historical comparisons can be affected by drift even with the same saved calibration.']
    lines += [
              '- Both raw and decoded errors are reported. Decoder-introduced errors can identify a mismatched DEM.',
              '- These known-state basis checks are not full state tomography and do not prove an arbitrary logical state was prepared correctly.']
    if comparisons:
        lines += ['', f'## {storage_us:g} us minus 15 us', '', '| Arm | Difference | 95% difference CI | Comparator |', '|---|---:|---|---|']
        for c in comparisons:
            lo, hi = c['difference_ci95']
            lines.append(f'| {c["source_arm"]} | {c["diagnostic_minus_15us"]:+.2%} | [{lo:+.2%}, {hi:+.2%}] | {c["provenance"]} |')
    (output/'analysis.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines)); print(f'Saved: {output}')


def execute(output, plan, timeout):
    if plan.get('diagnostic') not in ('zero_round_near_zero_storage_v2', 'zero_round_fixed_storage_v1'):
        raise RuntimeError('Wrong diagnostic plan; prepare the two-arm version')
    cases = plan['arms']; meta = plan['metadata']; layout = {int(q): v for q, v in meta['layout'].items()}
    validate(cases, layout)
    assert all(a.get('target_storage_ns', 0.) == meta['target_post_preparation_ns'] for a in cases)
    rows = memory.compile_durations(memory.get_compiler(meta), cases)
    (output/'pre_submission_timing.json').write_text(json.dumps(rows, indent=2))
    verify_timing(cases, rows, meta['max_near_zero_ns'])
    for a in cases:
        assert signature(a['prep_ops']) == a['preparation_sha256'] and signature(a['readout_ops']) == a['readout_sha256']
    from iqm.iqm_client import IQMClient
    client = IQMClient('https://resonance.iqm.tech', quantum_computer='emerald')
    dqa = client.get_dynamic_quantum_architecture(UUID(meta['calibration_set_id']))
    memory.check_mapping(cases, dqa.model_dump(mode='json'))
    for block in range(meta['blocks']):
        ordered = cases[block % len(cases):] + cases[:block % len(cases)]
        paths = [output/f'block{block}_{a["name"]}_raw.json' for a in ordered]
        if all(p.exists() for p in paths): raws = [json.loads(p.read_text()) for p in paths]
        else:
            job_path = output/f'block{block}_job.json'
            if job_path.exists():
                saved = json.loads(job_path.read_text()); assert saved['arms'] == [a['name'] for a in ordered]
                job = client.get_job(UUID(saved['job_id']))
            else:
                print(f'Submitting block {block}: {len(ordered)} circuits x {meta["shots_per_arm_per_block"]} shots', flush=True)
                job = client.submit_circuits(memory.native_circuits(ordered), shots=meta['shots_per_arm_per_block'], calibration_set_id=UUID(meta['calibration_set_id']))
                job_path.write_text(json.dumps(dict(job_id=str(job.job_id), arms=[a['name'] for a in ordered]), indent=2))
            raws = memory.collect_job_results(job, len(ordered), timeout, output/f'block{block}_status.json')
            for p, raw in zip(paths, raws): p.write_text(json.dumps(raw))
        for a, raw in zip(ordered, raws):
            directory = output/a['name']; ref = stim.Circuit.from_file(str(directory/'reference.stim'))
            measurements, dets, obs = memory.extract(raw, ref)
            assert len(obs) == meta['shots_per_arm_per_block']
            if a['kind'] == 'single': pred = np.zeros_like(obs)
            else:
                dem = stim.DetectorErrorModel.from_file(str(directory/'model.dem'))
                matcher = pymatching.Matching.from_detector_error_model(dem, enable_correlations=meta['correlated'])
                pred = matcher.decode_batch(dets, enable_correlations=meta['correlated'])
            np.savez_compressed(directory/f'block{block}.npz', raw_meas=measurements, det_events=dets, obs_flips=obs, predictions=pred)
    analyze(output)


def main(storage_ns=0., results_folder='zero_round_preparation'):
    p = argparse.ArgumentParser(description=__doc__ if storage_ns == 0 else f'Encoded logical X/Z, zero rounds, {storage_ns/1000:g} us storage; complete original encoder and readout.')
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--self-test', action='store_true'); mode.add_argument('--prepare', action='store_true')
    mode.add_argument('--hardware', action='store_true')
    p.add_argument('--source-run', type=Path, default=SOURCE); p.add_argument('--resume-dir', type=Path)
    p.add_argument('--shots', type=int, default=1000); p.add_argument('--blocks', type=int, default=2)
    p.add_argument('--max-near-zero-ns', type=float, default=8.); p.add_argument('--timeout-secs', type=float, default=600.)
    args = p.parse_args()
    if args.shots < 1 or args.blocks < 1 or not 0 <= args.max_near_zero_ns <= 100 or args.timeout_secs <= 0: p.error('Invalid shots, blocks, near-zero limit or timeout')
    if args.hardware:
        if not args.resume_dir: p.error('--hardware requires --resume-dir')
        plan = json.loads((args.resume_dir/'plan.json').read_text())
        if plan['metadata'].get('target_post_preparation_ns', 0.) != storage_ns:
            p.error(f'This entry point requires a {storage_ns:g} ns storage plan; use the matching prepared directory')
        execute(args.resume_dir, plan, args.timeout_secs); return
    if args.resume_dir: p.error('--resume-dir is only for hardware')
    source = json.loads((args.source_run/'plan.json').read_text()); meta = deepcopy(source['metadata'])
    if meta['device'] != 'emerald': p.error('Expected an Emerald source plan')
    cases = make_cases(source, storage_ns)
    layout = {int(q): v for q, v in meta['layout'].items()}; validate(cases, layout)
    if args.self_test:
        print(f'PASS: source encoder/readout unchanged; {storage_ns:g} ns storage; zero intermediate measurements.'); return
    meta.update(source_run=str(args.source_run.resolve()), shots_per_arm_per_block=args.shots, blocks=args.blocks,
                target_post_preparation_ns=float(storage_ns), max_near_zero_ns=args.max_near_zero_ns,
                with_15us_controls=False, diagnostic=f'Encoded logical X/Z only, zero extraction rounds, {storage_ns:g} ns storage')
    old_results = args.source_run/'results.json'
    if old_results.exists():
        previous = json.loads(old_results.read_text())['summary']
        meta['historical_15us_summary'] = {name: previous[name] for name in ARM_NAMES}
    rows = memory.compile_durations(memory.get_compiler(meta), cases)
    verify_timing(cases, rows, args.max_near_zero_ns)
    output = HERE.parent/'output'/results_folder/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    output.mkdir(parents=True)
    plan = dict(diagnostic='zero_round_near_zero_storage_v2' if storage_ns == 0 else 'zero_round_fixed_storage_v1', metadata=meta, arms=cases, timing=rows)
    (output/'plan.json').write_text(json.dumps(plan, indent=2))
    for a in cases:
        directory = output/a['name']; directory.mkdir()
        memory.model(a).to_file(str(directory/'reference.stim'))
        noisy = memory.model(a, True, meta['idle_error_per_us']); noisy.to_file(str(directory/'model.stim'))
        noisy.detector_error_model(decompose_errors=True).to_file(str(directory/'model.dem'))
        (directory/'native_instructions.json').write_text(json.dumps(memory.instructions(a), indent=2))
    lines = [f'# {storage_ns/1000:g} us storage: preparation and readout retained', '',
             '| Arm | Storage (ns) | Initialization (ns) | Encoder (ns) | Readout (ns) | CZ count |', '|---|---:|---:|---:|---:|---:|']
    for row in rows:
        lines.append(f'| {row["name"]} | {row["storage_before_readout_ns"]:.1f} | {row["initialization_ns"]:.1f} | {row["state_preparation_ns"]:.1f} | {row["final_readout_ns"]:.1f} | {row.get("compiled_cz_count", "unknown")} |')
    lines += ['', 'No hardware submitted. Initialization, encoder and final readout remain outside the storage window.']
    (output/'timing.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines)); print(f'Prepared {len(cases)} arms x {args.blocks} blocks x {args.shots} shots. Saved: {output}')


if __name__ == '__main__': main()
