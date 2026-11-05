"""Exercise interrupted refinement installation with clearly marked fake data.

This checks file integrity/resumption only. It never starts COMSOL or claims a
synthetic fixture is a solved FE model.
"""
import hashlib,json,time
from pathlib import Path
import run_transient as rt
from refine_bulk_transient import prepare_case,commit_trial,PAIRS,read


def main():
    root=rt.ROOT/'output/transient_refined/bulk_refinement_checks'/f'attempt_{time.time_ns()}'
    root.mkdir(parents=True);(root/'FIXTURE_ONLY.txt').write_text('Synthetic bookkeeping fixtures, not simulation results.')
    for variant,case in PAIRS:
        p=root/'bulk'/variant/f'L8/seed_101/{case}';(p/'blocks/P').mkdir(parents=True)
        signature=dict(layout=dict(kind='bulk',profile='refined_v6',variant=variant,dt=1. if variant=='base' else .5,end_time=10.),frames=case=='group')
        record=dict(signature=signature,status='complete',fingerprint=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest(),runs=[])
        rt.write_json(p/'run.json',record);(p/'P_history.npz').write_bytes(b'fixture: original history, not a real npz')
        (p/'P_solve.json').write_text('{}');(p/'blocks/P/restart_002.mph').write_bytes(b'fixture: opaque checkpoint, not a real mph')
        rt.write_json(p/'blocks/P/progress.json',dict(completed=3,checkpoint_file='restart_002.mph',window_converged=True,stable_windows=2))
    # Check promoted/extended signatures against fresh physical case requests,
    # including the actual reference-window policy, without solving fake data.
    cfg=read(rt.ROOT/'config.json');C,rho=rt.validate_material(cfg)
    def reference_layout(variant):
        s=rt.layout('bulk',8.,5e6,C,rho,variant,'refined_v6')
        arrival=(s['positions'][-1]-s['source_x'])/(C[0,0]/rho)**.5+1e-6
        s.update(end_time=max(50e-6,4*arrival),adaptive_window=True,minimum_window=12e-6)
        return s
    probe=root/'signature_probe';probe.mkdir()
    signature=dict(schema=rt.SCHEMA,layout=reference_layout('time'),material=cfg,case='reference',seed=101,
                   geometry_hash=hashlib.sha256(b'').hexdigest(),frames=False)
    rt.write_json(probe/'run.json',dict(status='complete',signature=signature,runs=[],
        fingerprint=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()))
    promoted_probe=root/'signature_promoted';prepare_case(probe,promoted_probe,'promote_time')
    expected=dict(signature,layout=reference_layout('base'))
    expected['layout'].update(dt=expected['layout']['dt']/2,time_step_refinement=2.)
    assert read(promoted_probe/'run.json')['signature']==expected
    (promoted_probe/'blocks/P').mkdir(parents=True)
    rt.write_json(promoted_probe/'blocks/P/progress.json',dict(completed=0))
    extended_probe=root/'signature_extended';prepare_case(promoted_probe,extended_probe,'extend_window')
    expected['layout'].update(end_time=expected['layout']['end_time']*2,window_extension_factor=2.)
    assert read(extended_probe/'run.json')['signature']==expected
    trial=root/'bulk_trials/attempt_001';source=root/'bulk/time/L8/seed_101/reference';promoted=trial/'bulk/base/L8/seed_101/reference'
    prepare_case(source,promoted,'promote_time')
    target=read(promoted/'run.json');original=read(source/'run.json')
    assert target['signature']['layout']['dt']==original['signature']['layout']['dt']
    assert target['signature']['layout']['variant']=='base' and target['signature']['layout']['time_step_refinement']==2
    assert (promoted/'P_history.npz').read_bytes()==(source/'P_history.npz').read_bytes()
    source=root/'bulk/base/L8/seed_101/group';extended=trial/'bulk/base/L8/seed_101/group'
    prepare_case(source,extended,'extend_window')
    assert not (extended/'P_history.npz').exists() and not (extended/'P_solve.json').exists()
    state=read(extended/'blocks/P/progress.json');assert state['completed']==3 and not state['window_converged'] and state['stable_windows']==0
    assert read(source/'blocks/P/progress.json')['window_converged']
    state['completed']=4;rt.write_json(extended/'blocks/P/progress.json',state)
    prepare_case(source,extended,'extend_window');assert read(extended/'blocks/P/progress.json')['completed']==4
    for v,c in PAIRS:
        p=trial/'bulk'/v/f'L8/seed_101/{c}'
        if not (p/'run.json').exists():prepare_case(root/'bulk'/v/f'L8/seed_101/{c}',p,'copy')
        record=read(p/'run.json');record['status']='complete';rt.write_json(p/'run.json',record)
    plan=dict(attempt=1,trial=str(trial),fingerprints={f'{v}/{c}':read(trial/'bulk'/v/f'L8/seed_101/{c}/run.json')['fingerprint'] for v,c in PAIRS})
    archive=root/'bulk_history/attempt_001';archive.parent.mkdir()
    assert (root/'bulk').resolve().is_relative_to(root.resolve()) and archive.resolve().is_relative_to(root.resolve())
    (root/'bulk').rename(archive)  # Simulate interruption after only the first rename.
    commit_trial(root,plan);commit_trial(root,plan)
    assert (archive/'base/L8/seed_101/group/P_history.npz').read_bytes()==b'fixture: original history, not a real npz'
    lineage=read(root/'bulk/base/L8/seed_101/reference/bulk_refinement_provenance.json')
    assert Path(lineage['archived_source']).exists() and read(Path(lineage['archived_source'])/'run.json')['fingerprint']==lineage['source_fingerprint']
    rt.write_json(root.parent/'bookkeeping_validation.json',dict(passed=True,fixture=str(root),
        checks=['role promotion retains solved data','new signatures match fresh physical requests','window extension preserves checkpoint and resets early-stop state',
                'resuming preparation does not overwrite new progress','directory swap recovers after interruption',
                'installation is idempotent','original histories remain archived'],scope='Bookkeeping only; not a numerical wave validation.'))
    print('Bulk refinement bookkeeping checks passed; no COMSOL solve was started.')


if __name__=='__main__':main()
