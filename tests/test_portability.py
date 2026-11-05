import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import detectability_200kHz_L48 as d
import detectability_followon as follow
from cfrp_pores.portability import ROOT, project_path, runtime_root, atomic_savez, validate_npz


class PortableTests(unittest.TestCase):
    def test_negative_single_response_uses_absolute_score(self):
        self.assertEqual(d.classify(-.002, .0001, single=True), 'above')
        self.assertEqual(d.classify(-.002, .0001), 'below')
        self.assertEqual(d.classify(-.0005, .0001, single=True), 'below')
        self.assertEqual(d.classify(-.001, .0001, single=True), 'indeterminate')
        for loss, uncertainty in [(float('nan'), .001), (.001, -.001), (.001, float('inf'))]:
            with self.assertRaises(ValueError):
                d.classify(loss, uncertainty, single=True)

    def test_old_absolute_paths_always_resolve_inside_checkout(self):
        manifest = json.loads((ROOT/'migration/source_snapshot.json').read_text(encoding='utf-8'))
        old = manifest['origin'].replace('\\', '/')+'/output/sample/geometry.json'
        self.assertEqual(project_path(old), ROOT/'output/sample/geometry.json')
        self.assertEqual(project_path(old.replace('/', '\\')), ROOT/'output/sample/geometry.json')
        with self.assertRaises(ValueError):
            project_path(manifest['origin']+'/../outside.json')

    def test_runtime_is_configurable_without_e_drive(self):
        with tempfile.TemporaryDirectory() as name, patch.dict('os.environ', {'CFRP_RUNTIME_DIR': name}):
            self.assertEqual(runtime_root(), Path(name).resolve())
        with patch.dict('os.environ', {'CFRP_RUNTIME_DIR': 'D:/含中文/runtime'}):
            with self.assertRaises(ValueError):
                runtime_root()

    def test_atomic_arrays_preserve_previous_cache_on_failed_write(self):
        with tempfile.TemporaryDirectory() as name:
            path = Path(name)/'arrays.npz'
            atomic_savez(path, value=np.arange(4))
            before = path.read_bytes()
            with self.assertRaises(ValueError):
                atomic_savez(path, value=np.array([float('nan')]))
            self.assertEqual(path.read_bytes(), before)
            self.assertEqual(list(Path(name).glob('*.tmp')), [])
            validate_npz(path, ('value',))
            with self.assertRaises(ValueError):
                validate_npz(path, ('missing',))

    def test_corrupt_array_and_modified_transient_are_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            folder = Path(name)
            (folder/'broken.npz').write_bytes(b'interrupted')
            with self.assertRaises(Exception):
                validate_npz(folder/'broken.npz')
            signature = dict(layout={})
            d.write(folder/'signature.json', signature)
            d.write(folder/'complete.json', dict(signature_sha256=d.signature(signature), history_sha256={'S0_history.npz':'invalid'}))
            (folder/'S0_history.npz').write_bytes(b'tampered')
            with self.assertRaisesRegex(RuntimeError, 'history changed'):
                follow.analyze_transient(folder)


if __name__ == '__main__':
    unittest.main()
