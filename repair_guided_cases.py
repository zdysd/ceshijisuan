"""Preserve invalid CAD outputs and rerun exactly the same pore contours."""
import mph,site
from pathlib import Path
ROOT=Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import shutil
from guided_collective import OUT,settings,grouped_pores,run_case,read

if __name__=='__main__':
    client=mph.Client(cores=4,version='6.2');s=settings()
    for d,seed in [(20,106),(40,105)]:
        folder=OUT/f'd{d}/L8/seed_{seed}'
        if read(folder/'results.json').get('signature',{}).get('absolute_repair_tolerance_m')==1e-9:continue
        backup=folder/'invalid_cad_before_repair';backup.mkdir(exist_ok=True)
        for name in ['model.mph','results.json','signals.npz']:
            if not (backup/name).exists():shutil.copy2(folder/name,backup/name)
        pores,_=grouped_pores(d,8,seed)
        run_case(client,folder,s,pores,force=True,repair_tolerance_m=1e-9)
