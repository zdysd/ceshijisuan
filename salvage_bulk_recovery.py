"""Commit a stopped P recovery prefix only after checking its native restart seam."""
import argparse
import hashlib
import json
import time

import numpy as np
import run_transient as rt


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--commit', action='store_true')
    args = parser.parse_args()
    folder = rt.ROOT/'output/transient_converged/bulk/base/L8/seed_101/group'
    blocks = folder/'blocks/P'
    progress_path = blocks/'progress.json'
    progress_hash = digest(progress_path)
    record = json.loads((folder/'run.json').read_text(encoding='utf-8-sig'))
    state = json.loads(progress_path.read_text(encoding='utf-8-sig'))
    s = record['signature']['layout']
    index = state['completed']
    checkpoint = blocks/state.get('checkpoint_file', 'restart.mph')
    checkpoint_hash = digest(checkpoint)
    previous = dict(np.load(blocks/f'history_{index-1:03d}.npz'))
    root = rt.ROOT/'output/transient_refined/recovery_salvage'/f'attempt_{time.time_ns()}'
    root.mkdir(parents=True)
    client = rt.start_session(root, 4)
    import jpype
    prior = client.load(str(checkpoint.resolve()))
    try:
        native_time = np.asarray(prior.evaluate('t')).ravel()
        np.testing.assert_allclose(native_time, previous['times'][-1:], rtol=1e-12, atol=1e-17)
        initial_u = np.asarray(prior.java.sol('sol1').getU(jpype.JInt(1))).copy()
        initial_v = np.asarray(prior.java.sol('sol1').getUDot(jpype.JInt(1))).copy()
    finally:
        client.remove(prior)
    model = rt.recover_mode(client, record, 'P')
    if model is None:
        raise RuntimeError('No stopped recovery model is available.')
    try:
        times = np.asarray(model.evaluate('t')).ravel()
        np.testing.assert_allclose(times[:1], previous['times'][-1:], rtol=1e-12, atol=1e-17)
        if len(times) < 2 or times[-1] <= times[0]:
            raise RuntimeError('Recovery has no new complete output interval.')
        np.testing.assert_allclose(np.diff(times), s['output_dt'], rtol=1e-10, atol=1e-17)
        np.testing.assert_allclose(times/s['output_dt'], np.round(times/s['output_dt']), rtol=0, atol=1e-8)
        for name in ['dt', 'tout', 'f0']:
            requested = s[{'dt':'dt', 'tout':'output_dt', 'f0':'f0'}[name]]
            np.testing.assert_allclose(model.java.param().evaluate(name), requested, rtol=1e-12, atol=0)
        sol = model.java.sol('sol1')
        errors = {}
        for name, getter, expected in [('displacement', sol.getU, initial_u),
                                       ('velocity', sol.getUDot, initial_v)]:
            actual = np.asarray(getter(jpype.JInt(1)))
            errors[name] = float(np.linalg.norm(actual-expected)/max(np.linalg.norm(expected), 1e-100))
        history = rt.extract_history(model, dict(s, end_time=float(times[-1])))
        C, rho = rt.validate_material(json.loads((rt.ROOT/'config.json').read_text()))
        scale = np.array([1, 1, 2*np.pi*s['f0']*np.sqrt(rho*C[0,0]),
                          2*np.pi*s['f0']*np.sqrt(rho*C[0,0])])[None,:,None]
        errors['receiver'] = float(np.linalg.norm((history['values'][0]-previous['values'][-1])/scale)
                                  /max(np.linalg.norm(previous['values'][-1]/scale), 1e-100))
        report = dict(source=record['active_recovery_directory'],
                      source_checkpoint_sha256=checkpoint_hash, old_completed=index,
                      start_time_s=float(times[0]), end_time_s=float(times[-1]),
                      output_times=len(times), seam_relative_errors=errors,
                      tolerance=1e-10, passed=max(errors.values()) < 1e-10, committed=False)
        rt.write_json(root/'validation.json', report)
        print(json.dumps(report), flush=True)
        if not report['passed']:
            raise RuntimeError('Recovered native displacement/velocity or receivers differ from the committed state.')
        if args.commit:
            # Retain history and real field frames before compacting native data.
            np.savez_compressed(root/'recovered_history.npz', **history)
            rt.extract_movie(model, s, root, 'P', times)
            final_u = np.asarray(sol.getU(jpype.JInt(len(times)))).copy()
            final_v = np.asarray(sol.getUDot(jpype.JInt(len(times)))).copy()
            staged = root/'restart.mph'
            compact = rt.compact_transient(model, staged, times)
            if not compact['derivative_retained']:
                raise RuntimeError('Recovery lost its final velocity.')
            client.remove(model)
            model = client.load(str(staged.resolve()))
            for getter, expected in [(model.java.sol('sol1').getU, final_u),
                                      (model.java.sol('sol1').getUDot, final_v)]:
                actual = np.asarray(getter(jpype.JInt(1)))
                np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-25)
            if digest(progress_path) != progress_hash or digest(checkpoint) != checkpoint_hash:
                raise RuntimeError('Committed source changed during recovery; refusing publication.')
            import shutil
            destination = blocks/f'restart_recovered_{index:03d}.mph'
            if destination.exists():
                raise RuntimeError('Recovered checkpoint already exists; inspect its prior attempt.')
            shutil.copy2(root/'recovered_history.npz', blocks/f'history_{index:03d}.npz')
            shutil.copy2(root/'P_fields.npz', blocks/f'P_{index:03d}_fields.npz')
            shutil.copy2(staged, destination)
            state['checks'].append(dict(block=index, seam_relative_error=errors['receiver'],
                                        velocity_relative_error=errors['velocity'],
                                        recovery_salvage=True))
            state.update(completed=index+1, checkpoint_file=destination.name)
            state.setdefault('recovery_salvages', []).append(str(root/'validation.json'))
            rt.write_json(progress_path, state)
            report.update(committed=True, checkpoint_file=str(destination), final_state_and_velocity_reloaded=True)
            rt.write_json(root/'validation.json', report)
            rt.write_json(root.parent/'latest.json', report)
            print('COMMITTED RECOVERY', json.dumps(report), flush=True)
    finally:
        if model is not None:
            client.remove(model)


if __name__ == '__main__':
    main()
