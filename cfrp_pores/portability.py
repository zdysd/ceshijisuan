"""Paths and durable arrays for a relocated Windows checkout; no solver startup."""
import hashlib
import json
import os
from pathlib import Path
import tempfile
import uuid

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def project_path(value):
    text = str(value).replace('\\', '/')
    manifest = ROOT / 'migration/source_snapshot.json'
    origin = json.loads(manifest.read_text(encoding='utf-8'))['origin'].replace('\\', '/')
    if text.casefold().startswith(origin.casefold().rstrip('/') + '/'):
        relative = text[len(origin.rstrip('/')) + 1:]
        target = (ROOT / relative).resolve()
        if not target.is_relative_to(ROOT):
            raise ValueError('Input path escapes repository')
        return target
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def runtime_root():
    default = Path(os.environ.get('LOCALAPPDATA', tempfile.gettempdir())) / 'cfrp_comsol_runtime'
    path = Path(os.environ.get('CFRP_RUNTIME_DIR', str(default))).resolve()
    if not path.is_absolute() or not str(path).isascii():
        raise ValueError('Set CFRP_RUNTIME_DIR to an absolute ASCII path on a disk with sufficient space')
    return path


def validate_npz(path, required=()):
    with np.load(path, allow_pickle=False) as data:
        if not data.files or not set(required).issubset(data.files):
            raise ValueError(f'Incomplete array cache: {path}')
        for key in data.files:
            value = data[key]
            if np.issubdtype(value.dtype, np.number) and not np.isfinite(value).all():
                raise ValueError(f'Nonfinite array {key}: {path}')


def atomic_savez(path, **arrays):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with temporary.open('wb') as handle:
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        validate_npz(temporary)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()
