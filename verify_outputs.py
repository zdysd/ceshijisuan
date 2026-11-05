"""Reopen solved COMSOL models, recompute evidence, and solve a matched single pore.

Run with the ASCII-path COMSOL Python runtime. Original models are read only;
new verification results and single-pore models go under output/revalidation.
"""
import mph
import site
from pathlib import Path
ROOT = Path(__file__).resolve().parent
site.addsitedir(str(ROOT / '.venv/Lib/site-packages'))
import csv
import json
import sys
import subprocess
import numpy as np
from cfrp_pores.comsol import (create_2d, add_load, add_frequency_study,
                              save_and_solve, integrate_section, interpolate)
from cfrp_pores.attenuation import fit_attenuation, forward_bulk_displacement
from cfrp_pores.dispersion import validate_material
from cfrp_pores.pores import load_polygons, simplify_geometry

sys.stdout.reconfigure(encoding='utf-8')


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def forward(model, positions, width, frequency, impedance):
    result = []
    for x in positions:
        values = integrate_section(model, float(x), width, ['u', 'solid.sx', '1'])
        if values[2].real <= 0:
            raise RuntimeError('Empty solid receiver section.')
        result.append(forward_bulk_displacement(values[0]/values[2].real,
                                               values[1]/width, frequency, impedance))
    result = np.asarray(result)
    if not np.isfinite(result).all():
        raise RuntimeError('Nonfinite projected solution.')
    return result


def project_file(client, record, impedance, downstream=False):
    path = Path(record['model'])
    model = client.load(str(path))
    try:
        if len(model.java.sol('sol1').getU()) != record['dofs']:
            raise RuntimeError(f'DOF count differs in {path}.')
        values = forward(model, record['positions_m'], record['width_m'],
                         record['frequency_hz'], impedance)
        post = None
        if downstream:
            post = forward(model, [.00145, .00155, .00165], record['width_m'],
                           record['frequency_hz'], impedance)
        return values, post
    finally:
        client.remove(model)


def main():
    out = ROOT/'output/revalidation'
    out.mkdir(parents=True, exist_ok=True)
    cfg = read(ROOT/'config.json')
    C, rho = validate_material(cfg)
    impedance = np.sqrt(C[0, 0]*rho)
    checks = dict(material_status=cfg['material_status'],
                  scope='Numerical verification; no experimental detectability threshold.')
    unit = subprocess.run([str(ROOT/'.venv/Scripts/python.exe'), '-m', 'unittest',
                           'discover', '-s', 'tests', '-v'], cwd=ROOT,
                          capture_output=True, text=True, encoding='utf-8', errors='replace')
    (out/'unit_checks.txt').write_text(unit.stdout+unit.stderr, encoding='utf-8')
    checks['unit_checks_passed'] = unit.returncode == 0
    if unit.returncode:
        raise RuntimeError('Scientific checks failed; see unit_checks.txt.')
    client = mph.Client(cores=2, version='6.2')
    ensembles = []
    arrays = {}
    contrasts = []
    for folder in ['final_ensemble_5MHz', 'final_ensemble_10MHz']:
        root = ROOT/'output'/folder
        reference = read(root/'reference/measurement.json')
        frequency = reference['frequency_hz']
        ref, ref_post = project_file(client, reference, impedance, downstream=True)
        signals, post_signals, seeds, individual_error = [], [], [], []
        with (root/'forward_individual.csv').open(encoding='utf-8-sig') as stream:
            previous = {int(r['seed']): float(r['alpha_np_m']) for r in csv.DictReader(stream)}
        for item in sorted(root.glob('seed_*')):
            record = read(item/'measurement.json')
            signal, post = project_file(client, record, impedance, downstream=True)
            fitted = fit_attenuation(record['positions_m'], signal, ref)
            signals.append(signal); post_signals.append(post); seeds.append(record['seed'])
            individual_error.append(abs(fitted['alpha_np_m']-previous[record['seed']]))
            print(f'REOPEN {folder} seed={record["seed"]}: alpha={fitted["alpha_np_m"]:.9g}', flush=True)
        signals = np.asarray(signals); post_signals = np.asarray(post_signals)
        positions = np.asarray(reference['positions_m'])
        fit = fit_attenuation(positions, signals.mean(axis=0), ref)
        old = read(root/'forward_summary.json')
        rng = np.random.default_rng(1234)
        boot = [fit_attenuation(positions, signals[rng.integers(len(seeds), size=len(seeds))].mean(axis=0),
                                ref)['alpha_np_m'] for _ in range(500)]
        ci = np.percentile(boot, [2.5, 97.5])
        ensembles.append(dict(frequency_hz=frequency, samples=len(seeds), fit=fit,
                              ci_95_np_m=ci.tolist(),
                              maximum_individual_alpha_reextraction_error_np_m=max(individual_error),
                              coherent_alpha_reextraction_error_np_m=abs(fit['alpha_np_m']-old['fit']['alpha_np_m']),
                              ci_reextraction_error_np_m=float(np.max(abs(ci-old['bootstrap_ci_95_np_m']))),
                              reproduces_saved_result=max(individual_error)<1e-6,
                              positive_local_coherent_attenuation=bool(ci[0]>0),
                              relative_ci_half_width=float((ci[1]-ci[0])/(2*fit['alpha_np_m']))))
        tag = str(int(frequency/1e6))
        arrays[f'signals_{tag}MHz'] = signals
        arrays[f'reference_{tag}MHz'] = ref
        arrays[f'positions_{tag}MHz_m'] = positions
        arrays[f'post_signals_{tag}MHz'] = post_signals
        arrays[f'post_reference_{tag}MHz'] = ref_post
        # Retain the largest pore from seed 101, at precisely its group location.
        first_record = read(root/'seed_101/measurement.json')
        polygons, _ = load_polygons(root/'seed_101/geometry/geometry.json')
        polygons, _ = simplify_geometry(polygons, first_record['contour_tolerance_um'])
        largest = max(polygons, key=lambda p: p.area)
        mesh = first_record['mesh_m']
        model = create_2d(client, f'single_pore_{tag}MHz', cfg, [largest],
                          size_m=(.002, .001), pore_offset_m=(.0003, .00004),
                          mesh_m=mesh, lateral_periodic=True)
        try:
            add_load(model, ['1e6[Pa]', '0', '0']); add_frequency_study(model, frequency)
            info = save_and_solve(model, out/f'single_pore_{tag}MHz.mph')
            single_post = forward(model, [.00145, .00155, .00165], .001, frequency, impedance)
        finally:
            client.remove(model)
        single_T = single_post/ref_post
        group_T = post_signals[0]/ref_post
        mean_T = post_signals.mean(axis=0)/ref_post
        single_change = float(np.mean(1-abs(single_T)))
        group_change = float(np.mean(1-abs(group_T)))
        contrasts.append(dict(frequency_hz=frequency, pore_count=len(polygons),
                              largest_pore_equivalent_diameter_um=float(np.sqrt(4*largest.area/np.pi)),
                              largest_pore_bbox_um=[float(largest.bounds[2]-largest.bounds[0]), float(largest.bounds[3]-largest.bounds[1])],
                              reference_wavelength_um=float(np.sqrt(C[0,0]/rho)/frequency*1e6),
                              single_pore_amplitude_reduction=single_change,
                              same_geometry_group_amplitude_reduction=group_change,
                              ensemble_group_amplitude_reduction=float(np.mean(1-abs(mean_T))),
                              group_response_larger=group_change>single_change,
                              single_pore_model=info,
                              observation='Forward transmission in the same pore-free downstream sections; no detectability threshold.',
                              boundary_assumption='Transverse periodic cell; single pore is also periodically repeated.',
                              warning='Different total porosities; this contrast measures adding pores, not a clustering or nonlinear-interaction effect.'))
        arrays[f'single_transmission_{tag}MHz'] = single_T
        arrays[f'group_transmission_{tag}MHz'] = group_T
        arrays[f'ensemble_transmission_{tag}MHz'] = mean_T
        print(f'CONTRAST {tag}MHz: single={single_change:.6g}, group={group_change:.6g}', flush=True)
    checks['ensembles'] = ensembles
    checks['single_group_contrasts'] = contrasts
    # Re-extract every saved discretization/domain comparison, without resaving.
    validation_rows = []
    for case in ['base', 'mesh_fine', 'contour_full', 'length_long', 'width_wide']:
        root = ROOT/'output/final_validation'/case
        reference = read(root/'reference/measurement.json'); pore = read(root/'pores/measurement.json')
        ref, _ = project_file(client, reference, impedance); signal, _ = project_file(client, pore, impedance)
        fit = fit_attenuation(pore['positions_m'], signal, ref)
        validation_rows.append(dict(case=case, **fit))
        print(f'DOMAIN {case}: alpha={fit["alpha_np_m"]:.9g}', flush=True)
    base = validation_rows[0]['alpha_np_m']
    for item in validation_rows:
        item['relative_alpha_change_from_base'] = abs(item['alpha_np_m']-base)/max(abs(base), abs(item['alpha_np_m']))
        item['relative_coefficient_passed_5percent'] = item['relative_alpha_change_from_base'] <= .05
    checks['coefficient_convergence'] = validation_rows
    checks['domain_coefficient_converged_5percent'] = all(r['relative_coefficient_passed_5percent'] for r in validation_rows)
    # Verify both exported COMSOL models actually contain the saved solved fields.
    root = ROOT/'output/equivalent_bulk_5MHz'
    model = client.load(str(root/'equivalent_P_5MHz.mph'))
    try:
        d = np.load(root/'comparison.npz')
        values = interpolate(model, [d['x_m']], ['u'])[0, -1]
        error = float(np.max(abs(values-d['comsol'])))
        checks['scalar_model_reopen'] = dict(max_saved_field_difference=error, passed=error<1e-10)
    finally:
        client.remove(model)
    root = ROOT/'output/comsol_3d_mesh_based'
    record = read(root/'results.json'); model = client.load(str(root/'random_pores.mph'))
    try:
        values = interpolate(model, [[.00012,.0003,.00048],[.00015]*3,[.00015]*3], ['u','v','w'])
        old = np.load(root/'receivers.npz')['displacements']
        error = float(np.max(abs(values-old))/max(float(np.max(abs(old))), 1e-30))
        dofs = len(model.java.sol('sol1').getU())
        elements = int(model.java.component('comp1').mesh('mesh1').getNumElem())
        checks['three_dimensional_reopen'] = dict(dofs=dofs, mesh_elements=elements,
                                                  relative_probe_difference=error,
                                                  passed=error<1e-9 and dofs==record['dofs'] and elements==record['mesh_elements'])
    finally:
        client.remove(model)
    checks['experimental_single_pore_undetectability_verified'] = False
    checks['centimetre_extrapolation_validated'] = False
    checks['numerical_reproduction_passed'] = bool(checks['unit_checks_passed'] and
        all(x['reproduces_saved_result'] for x in ensembles) and
        checks['scalar_model_reopen']['passed'] and checks['three_dimensional_reopen']['passed'])
    np.savez_compressed(out/'reextracted_signals.npz', **arrays)
    (out/'verification.json').write_text(json.dumps(checks, indent=2), encoding='utf-8')
    print(json.dumps(checks, indent=2), flush=True)


if __name__ == '__main__':
    main()
