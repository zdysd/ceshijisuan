"""Check CAD cavity counts and area on completed guided-study models."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import json,sys
from guided_collective import OUT,read
sys.stdout.reconfigure(encoding='utf-8')


def main():
    client=None;rows=[]
    audit=OUT/'geometry_audit.json'
    previous={r['case']:r for r in read(audit)} if audit.exists() else {}
    for path in sorted(OUT.rglob('model.mph')):
        if 'invalid_cad_before_repair' in path.parts:continue
        if not (path.parent/'results.json').exists():continue
        meta_path=path.parent/'geometry/metadata.json'
        expected=read(meta_path)['pore_count'] if meta_path.exists() else (1 if 'single' in path.parent.name or 'maximum' in path.parent.name else None)
        if meta_path.exists() and ('single' in path.parent.name):expected=1
        if expected is None:continue
        case=str(path.parent.relative_to(OUT));record=read(path.parent/'results.json')
        if 'finite_voids' in record:
            count=record['finite_voids'];source='Native final-geometry count recorded during this solve'
        elif case in previous and previous[case]['passed']:
            count=previous[case]['cad_cavities'];source='Previously inspected unchanged saved model'
        else:
            if client is None:client=mph.Client(cores=2,version='6.2')
            model=client.load(str(path))
            try:count=int(model.java.component('comp1').geom('geom1').getNFiniteVoids())
            finally:client.remove(model)
            source='Saved model inspected in this audit'
        row=dict(case=case,expected_cavities=expected,cad_cavities=count,passed=count==expected,source=source)
        rows.append(row);print(json.dumps(row),flush=True)
    (OUT/'geometry_audit.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')


if __name__=='__main__':main()
