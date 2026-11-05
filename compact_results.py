"""Compact only this project's generated transient models, preserving probes."""
import mph
import site
from pathlib import Path
site.addsitedir(str(Path(__file__).parent/'.venv/Lib/site-packages'))
import argparse
import json
import sys
import numpy as np
from cfrp_pores.comsol import compact_transient,interpolate
sys.stdout.reconfigure(encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('folder');args=parser.parse_args()
    workspace=Path(__file__).parent.resolve();folder=Path(args.folder).resolve()
    if not folder.is_relative_to(workspace/'output'):raise ValueError('Only generated models under this project output/ can be compacted.')
    client=mph.Client(cores=2,version='6.2');records=[]
    for history in folder.glob('*_transient.npz'):
        source=folder/(history.name.removesuffix('_transient.npz')+'.mph')
        if not source.exists():continue
        data=np.load(history)
        before=source.stat().st_size
        if before<20_000_000:continue
        model=client.load(str(source))
        temporary=source.with_name(source.stem+'_compact_verified.mph')
        info=compact_transient(model,temporary,data['times_s'])
        positions=data['receivers_m'];width=float(model.java.component('comp1').geom('geom1').getBoundingBox()[3])
        actual=interpolate(model,[positions,np.full(len(positions),width/2)])[:,0,:].real
        if not np.allclose(actual,data['displacement'][:,-1,:],rtol=1e-8,atol=1e-15):
            raise RuntimeError('Final state mismatch; the full model has not been replaced.')
        client.remove(model)
        reopened=client.load(str(temporary))
        actual=interpolate(reopened,[positions,np.full(len(positions),width/2)])[:,0,:].real
        if not np.allclose(actual,data['displacement'][:,-1,:],rtol=1e-8,atol=1e-15):
            raise RuntimeError('Reopened state mismatch; the full model has not been replaced.')
        client.remove(reopened)
        temporary.replace(source)
        record=dict(model=str(source),bytes_before=before,bytes_after=source.stat().st_size,
                    full_receiver_history_preserved=True,final_state_reopened_and_verified=True,**info)
        records.append(record);print(json.dumps(record),flush=True)
        (folder/'storage_compaction.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
