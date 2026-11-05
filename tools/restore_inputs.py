"""Restore the versioned compact snapshot without overwriting computed results."""
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]


def restore():
    manifest = json.loads((ROOT/'migration/portable_data_manifest.json').read_text(encoding='utf-8'))
    restored = 0
    for entry in manifest['files']:
        relative = Path(entry['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Unsafe snapshot path')
        source = ROOT/'data/portable_inputs'/relative
        if hashlib.sha256(source.read_bytes()).hexdigest() != entry['sha256']:
            raise ValueError(f'Portable input changed: {relative}')
        target = ROOT/relative
        if not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            restored += 1
    print(f'Verified {len(manifest["files"])} bundled files; restored {restored} missing files; existing results preserved.')


if __name__ == '__main__':
    restore()
