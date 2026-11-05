"""Export independent, solution-free COMSOL copies; never overwrite originals."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import argparse,json,hashlib,shutil,sys,time
import numpy as np
import jpype
from guided_collective import settings,read

DEST=ROOT/'output/COMSOL_reproduction_20261014'
SCHEMA='solution_free_selectors_v1'
sys.stdout.reconfigure(encoding='utf-8')


def inventory():
    rows=[]
    for p in sorted((ROOT/'output').rglob('*.mph')):
        if DEST in p.parents:continue
        rel=p.relative_to(ROOT/'output')
        main=rel.parts[0]=='guided_collective' and 'invalid_cad_before_repair' not in rel.parts
        category='main_models' if main else ('historical_invalid' if 'invalid_cad_before_repair' in rel.parts else 'historical_models')
        suffix=Path(*rel.parts[1:]) if main else rel
        rows.append(dict(source=str(p),relative_source=rel.as_posix(),category=category,
                         destination=(Path(category)/suffix).as_posix(),source_bytes=p.stat().st_size,source_mtime_ns=p.stat().st_mtime_ns))
    return rows


def choose(index,expressions):
    expression=expressions[-1]
    for i in reversed(range(len(expressions)-1)):
        expression=f'if({index}=={i+1},{expressions[i]},{expression})'
    return expression


def add_selectors(model,record,s):
    c=model.java.component('comp1');param=model.java.param()
    frequencies=record['signature']['frequencies']
    original_frequency=float(param.evaluate('f0'))
    original={n:np.array([[float(str(v)) for v in row] for row in c.func(n).getStringMatrix('table')]) for n in ['ur','ui','vr','vi']}
    default_index=frequencies.index(original_frequency)+1
    param.set('freq_id',str(default_index),'频率编号 见本模型允许的频率列表')
    param.set('mode_id','1','入射模态 0为S0 1为A0 其他值不受支持')
    param.set('u_amp','1e-8[m]','源位移幅值 对应原始仿真10纳米')
    param.set('f0',choose('freq_id',[f'{f:.16g}[Hz]' for f in frequencies]),
              '由freq_id选取 '+', '.join(f'{i+1}={f/1000:g}kHz' for i,f in enumerate(frequencies)))
    functions={};default_error=0.
    for mode_id,kind in enumerate(['S0','A0']):
        for fi,freq in enumerate(frequencies,1):
            wave=s['plate'].at_frequency(freq,kind)
            profiles={'ur':np.conj(wave['u']).real,'ui':np.conj(wave['u']).imag,
                      'vr':np.conj(wave['v']).real,'vi':np.conj(wave['v']).imag}
            for part,values in profiles.items():
                tag=f'pkg{mode_id}f{fi}{part}'
                fn=c.func().create(tag,'Interpolation');fn.label(f'{kind} {freq/1000:g}kHz {part}')
                fn.set('source','table');fn.set('funcname',tag)
                table=np.column_stack([wave['y_m']+s['h']/2,values])
                fn.set('table',jpype.JArray(jpype.JString,2)([[f'{a:.17g}',f'{b:.17g}'] for a,b in table]))
                fn.set('argunit','m');fn.set('fununit','1')
                functions[(mode_id,fi,part)]=tag+'(y)'
                if kind=='A0' and freq==original_frequency:
                    default_error=max(default_error,float(np.max(abs(table-original[part]))))
    components=[]
    for real,imag in [('ur','ui'),('vr','vi')]:
        modes=[]
        for mi in [0,1]:
            values=[f'({functions[mi,fi,real]}+i*{functions[mi,fi,imag]})' for fi in range(1,len(frequencies)+1)]
            modes.append(choose('freq_id',values))
        components.append('u_amp*if(mode_id==0,'+modes[0]+','+modes[1]+')')
    c.physics('solid').feature('source').set('U0',components+['0'])
    # Original interpolants remain embedded as an audit reference, but the source
    # now uses the explicit selectors and correctly matched modal shapes.
    if default_error>1e-6:raise ValueError(f'New default mode differs from source: {default_error}')
    return dict(frequencies_hz=frequencies,default_freq_id=default_index,default_mode_id=1,
                original_source_profile_max_difference=default_error,source_expressions=components)


def parameter_snapshot(model):
    p=model.java.param();out={}
    for name in p.varnames():
        name=str(name)
        out[name]=dict(expression=str(p.get(name)),description=str(p.descr(name)))
    return out


def clear_solutions(model):
    j=model.java
    for tag in j.sol().tags():j.sol(tag).clearSolutionData()
    for tag in j.result().table().tags():j.result().table(tag).clearTableData()
    # Numerical interpolation nodes can cache evaluations; retain dataset and
    # plot definitions but rebuild numerical evaluations after the new solve.
    for tag in list(j.result().numerical().tags()):j.result().numerical().remove(tag)
    return {str(tag):bool(j.sol(tag).isEmpty()) for tag in j.sol().tags()}


def copy_external_inputs(model,target):
    dependencies=[]
    def visit(feature,label):
        try:properties=feature.properties()
        except AttributeError:properties=[]
        for key in properties:
            key=str(key)
            if 'filename' not in key.lower() or key.lower().endswith('hidden'):continue
            try:value=str(feature.getString(key))
            except Exception:continue
            if not value:continue
            source=Path(value)
            if not source.is_absolute():source=ROOT/source
            if source.is_file():
                folder=target.parent/'assets';folder.mkdir(exist_ok=True)
                name=hashlib.sha256(str(source).encode()).hexdigest()[:10]+'_'+source.name
                copied=folder/name
                if not copied.exists():shutil.copy2(source,copied)
                feature.set(key,'assets/'+name)
                dependencies.append(dict(node=label,property=key,original=str(source),relative_path='assets/'+name,exists=True))
            else:dependencies.append(dict(node=label,property=key,original=value,exists=False))
        try:children=list(feature.feature().tags())
        except Exception:children=[]
        for tag in children:visit(feature.feature(tag),label+'/'+str(tag))
    j=model.java
    for ci in j.component().tags():
        c=j.component(ci)
        for category in ['geom','mesh','func']:
            items=getattr(c,category)()
            for tag in items.tags():visit(getattr(c,category)(tag),str(ci)+'/'+category+'/'+str(tag))
    for tag in j.func().tags():visit(j.func(tag),'func/'+str(tag))
    if dependencies:j.modelPath(str(target.parent))
    return dependencies


def export_one(client,row,s):
    source=Path(row['source']);target=DEST/row['destination'];info=target.with_suffix('.parameters.json')
    if target.exists() and info.exists():
        old=read(info)
        if old.get('schema')==SCHEMA and old.get('source_mtime_ns')==row['source_mtime_ns'] and old.get('reopened_solution_empty'):
            return old
    target.parent.mkdir(parents=True,exist_ok=True)
    print('LOAD',row['relative_source'],flush=True)
    model=client.load(str(source));extra={}
    try:
        if row['category']=='main_models':
            record=read(source.parent/'results.json')
            extra=add_selectors(model,record,s)
            extra.update(original_signature=record['signature'],material_status=record['material_status'],
                         finite_voids=record.get('finite_voids'))
            if (source.parent/'geometry/metadata.json').exists():extra['pore_metadata']=read(source.parent/'geometry/metadata.json')
        else:extra['external_inputs']=copy_external_inputs(model,target)
        parameters=parameter_snapshot(model)
        mesh={str(c):{str(m):int(model.java.component(c).mesh(m).getNumElem()) for m in model.java.component(c).mesh().tags()} for c in model.java.component().tags()}
        structure={str(c):dict(geometry_dimensions={str(g):int(model.java.component(c).geom(g).getSDim()) for g in model.java.component(c).geom().tags()},
                              physics={str(p):str(model.java.component(c).physics(p).getType()) for p in model.java.component(c).physics().tags()}) for c in model.java.component().tags()}
        states=clear_solutions(model)
        if not all(states.values()):raise RuntimeError('Solution data was not cleared.')
        model.java.label(target.name)
        model.save(str(target))
    finally:client.remove(model)
    check=client.load(str(target))
    try:
        empty={str(tag):bool(check.java.sol(tag).isEmpty()) for tag in check.java.sol().tags()}
        if not all(empty.values()):raise RuntimeError('Reopened model contains a solution.')
        assert parameter_snapshot(check)==parameters
        if row['category']=='main_models':
            for i,freq in enumerate(extra['frequencies_hz'],1):
                check.java.param().set('freq_id',str(i))
                assert abs(check.java.param().evaluate('f0')-freq)<1e-7
    finally:client.remove(check)
    result=dict(row,schema=SCHEMA,parameters=parameters,mesh_elements=mesh,structure=structure,
                cleared_solutions=states,reopened_solution_empty=all(empty.values()),
                exported_bytes=target.stat().st_size,**extra)
    info.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print('SAVED',row['destination'],target.stat().st_size,flush=True)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--scope',choices=['main','all','history'],default='main');p.add_argument('--limit',type=int);args=p.parse_args()
    DEST.mkdir(parents=True,exist_ok=True)
    items=inventory();(DEST/'source_inventory.json').write_text(json.dumps(items,ensure_ascii=False,indent=2),encoding='utf-8')
    if args.scope=='main':items=[r for r in items if r['category']=='main_models']
    if args.scope=='history':items=[r for r in items if r['category']!='main_models']
    if args.limit:items=items[:args.limit]
    client=mph.Client(cores=2,version='6.2');s=settings();rows=[]
    for index,row in enumerate(items,1):
        result=export_one(client,row,s);rows.append(result)
        (DEST/f'export_{args.scope}.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
        print('PROGRESS',index,len(items),flush=True)


if __name__=='__main__':main()
