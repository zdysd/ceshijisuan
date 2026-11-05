"""Offline preflight. Importing mph does not start a COMSOL client or solver."""
import argparse
import importlib
import json
import os
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comsol-home', type=Path, help='Check an installed COMSOL directory without launching it')
    args = parser.parse_args()
    checks = []
    checks.append(dict(name='python_3_12', passed=sys.version_info[:2] == (3, 12), value=sys.version))
    checks.append(dict(name='windows_solver_worker', passed=os.name == 'nt', value=os.name))
    for name in ['numpy','scipy','matplotlib','skimage','shapely','mph','jpype','meshpy']:
        try:
            module = importlib.import_module(name)
            checks.append(dict(name=name, passed=True, value=getattr(module, '__version__', 'imported')))
        except Exception as error:
            checks.append(dict(name=name, passed=False, value=str(error)))
    try:
        from cfrp_pores.portability import runtime_root
        path = runtime_root()
        path.mkdir(parents=True, exist_ok=True)
        import tempfile
        with tempfile.TemporaryFile(dir=path):
            pass
        checks.append(dict(name='runtime_path_writable_ascii', passed=True, value=str(path), free_gib=round(shutil.disk_usage(path).free/2**30, 2)))
    except Exception as error:
        checks.append(dict(name='runtime_path_writable_ascii', passed=False, value=str(error)))
    try:
        import detectability_worker as worker
        worker.verify_inputs()
        worker.base_settings()
        checks.append(dict(name='frozen_inputs_and_engine', passed=True))
    except Exception as error:
        checks.append(dict(name='frozen_inputs_and_engine', passed=False, value=str(error)))
    if args.comsol_home:
        checks.append(dict(name='comsol_installation_path', passed=args.comsol_home.is_dir(), value=str(args.comsol_home)))
    result = dict(offline_passed=all(v['passed'] for v in checks), checks=checks,
                  native_solver_tested=False, license_tested=False,
                  next_step='On the target computer only: verify COMSOL 6.2 build, licensed modules, then explicitly authorize a native smoke solve.')
    print(json.dumps(result, ensure_ascii=True, indent=2))
    return 0 if result['offline_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
