"""Incident-normalized transmission: pore-free sections before/after the ROI."""
import mph
import site
from pathlib import Path
ROOT = Path(__file__).resolve().parent
site.addsitedir(str(ROOT/'.venv/Lib/site-packages'))
import json
import sys
import numpy as np
from verify_outputs import read, forward
from cfrp_pores.comsol import create_2d, add_load, add_frequency_study, save_and_solve
from cfrp_pores.dispersion import validate_material
from cfrp_pores.pores import load_polygons, simplify_geometry
sys.stdout.reconfigure(encoding='utf-8')


def section_signals(client, path, f, impedance):
    model = client.load(str(path))
    try:
        return forward(model, [.00015, .00145, .00155, .00165], .001, f, impedance)
    finally:
        client.remove(model)


def main():
    out = ROOT/'output/revalidation'
    cfg = read(ROOT/'config.json'); C, rho = validate_material(cfg)
    impedance = np.sqrt(rho*C[0,0]); client = mph.Client(cores=2, version='6.2')
    results, arrays = [], {}
    for tag in ['5', '10']:
        f = float(tag)*1e6; root = ROOT/f'output/final_ensemble_{tag}MHz'
        ref = section_signals(client, root/'reference/model.mph', f, impedance)
        normalized_reference = ref[1:]/ref[0]
        group = []
        for folder in sorted(root.glob('seed_*')):
            values = section_signals(client, folder/'model.mph', f, impedance)
            group.append((values[1:]/values[0])/normalized_reference)
            print('INCIDENT NORMALIZATION', tag, folder.name, flush=True)
        group = np.asarray(group)
        polygons, _ = load_polygons(root/'seed_101/geometry/geometry.json')
        record = read(root/'seed_101/measurement.json')
        polygons, _ = simplify_geometry(polygons, record['contour_tolerance_um'])
        sorted_pores = sorted(polygons, key=lambda p:p.area)
        median = sorted_pores[len(sorted_pores)//2]
        model = create_2d(client, f'median_single_{tag}MHz', cfg, [median],
                          size_m=(.002,.001), pore_offset_m=(.0003,.00004),
                          mesh_m=record['mesh_m'], lateral_periodic=True)
        try:
            add_load(model, ['1e6[Pa]','0','0']); add_frequency_study(model,f)
            median_info = save_and_solve(model, out/f'median_single_{tag}MHz.mph')
        finally:
            client.remove(model)
        row = dict(frequency_hz=f, largest_pore_diameter_um=float(np.sqrt(4*sorted_pores[-1].area/np.pi)),
                   median_pore_diameter_um=float(np.sqrt(4*median.area/np.pi)),
                   group_pore_count=len(polygons), median_single_model=median_info,
                   observation='Cross-section coherent forward displacement; normalize by each model incident forward wave.',
                   receivers_m=[.00015,.00145,.00155,.00165],
                   periodic_cell=True, energy_transmittance=False,
                   experimental_detectability_verified=False)
        for label, path in [('largest_single',out/f'single_pore_{tag}MHz.mph'),
                            ('median_single',out/f'median_single_{tag}MHz.mph')]:
            values = section_signals(client,path,f,impedance)
            transmission = (values[1:]/values[0])/normalized_reference
            row[label+'_amplitude_reduction'] = float(np.mean(1-abs(transmission)))
            row[label+'_phase_shift_rad'] = float(np.angle(transmission.mean()))
            arrays[f'{label}_T_{tag}MHz'] = transmission
        row['group_seed101_amplitude_reduction'] = float(np.mean(1-abs(group[0])))
        mean = group.mean(axis=0)
        row['coherent_group_amplitude_reduction'] = float(np.mean(1-abs(mean)))
        row['coherent_group_phase_shift_rad'] = float(np.angle(mean.mean()))
        row['post_section_variation'] = float(np.ptp(abs(mean)))
        rng = np.random.default_rng(2026)
        boot = [float(np.mean(1-abs(group[rng.integers(len(group),size=len(group))].mean(axis=0)))) for _ in range(1000)]
        row['amplitude_reduction_ci95'] = np.percentile(boot,[2.5,97.5]).tolist()
        row['group_change_exceeds_median_single'] = row['coherent_group_amplitude_reduction']>row['median_single_amplitude_reduction']
        row['group_change_exceeds_largest_single'] = row['coherent_group_amplitude_reduction']>row['largest_single_amplitude_reduction']
        arrays[f'group_T_{tag}MHz'] = group; arrays[f'reference_{tag}MHz'] = ref
        results.append(row)
        print(json.dumps(row,indent=2),flush=True)
    result = dict(results=results, caveat='Transmission contrast is a collective response; it does not establish nonlinear interaction, experimental single-pore invisibility, or a material attenuation coefficient.')
    (out/'incident_normalized_transmission.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    np.savez_compressed(out/'incident_normalized_transmission.npz',**arrays)


if __name__=='__main__':
    main()
