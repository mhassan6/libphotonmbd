
import sys
sys.path.append("../src")
import numpy as np
import json
from atomic_constants import P_Ang

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode

from pytest import approx

pt = photon_mode(1, 0.2, 0.1, [0,0,1])
#filename = 'json/ar_dummy_periodic_unitcell.json'
filename = 'json/dummy_Ar_natm_2_dist_4_2k.json'
#filename = 'json/dummy_Ar_natm_1_dist_4_4k.json'

gauge = 'length'
output_array = pypmbd_parse_json(filename, pt, gauge)

vdwene_len = output_array['ene']
E0 = output_array['E0']
AT_energy = output_array['AT_energy']
rpa_orders = output_array['rpa_orders']
pt_AT_energy = output_array['pt_AT_energy']
pt_rpa_orders = output_array['pt_rpa_orders']
pt_number = output_array['pt_number']
pt_adadaa = output_array['pt_adadaa']

print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)

print('ene', vdwene_len, 'AT ene', AT_energy, 'ratio', vdwene_len/AT_energy)
print('log AT ene', AT_energy, 'rpa_order', np.sum(rpa_orders))

print('\n')
print('with photons')
print('ene', pt.exc, 'AT ene', pt_AT_energy, 'ratio', pt.exc/pt_AT_energy)
print('log AT ene', pt_AT_energy, 'rpa_order', np.sum(pt_rpa_orders))
print('photonnumber', pt_number)
print('photon adadaa', pt_adadaa)

#Compare first to MBD

from pymbd.fortran import MBDGeom

with open(filename, "r") as f:
  data = json.load(f)
  atomic_positions = P_Ang*np.array(data['atomic_positions']) #atomic_positions in a.u. while json in AA
  lattice = data['lattice'] #FIX to a.u.
  lattice = P_Ang*np.array(lattice) #lattice in a.u. while json in AA
  kpoints = data['kpoints']

print(lattice)
print(kpoints)

geom = MBDGeom(atomic_positions, lattice, kpoints, do_rpa=True, get_rpa_orders=True)
mbd_ene, mbd_orders = geom.mbd_energy_species(['Ar', 'Ar'], [1, 1], 0.83)

print('mbd_ene', mbd_ene)


print()
print('here')
print()

from pymbd import ang, from_volumes

def argon_crystal():
    return (
        np.array([(0.3, 0.1, 0.2), (4.1, -0.2, -0.1)]) * ang,
        np.array([(8.1, 0.1, -0.2), (0.3, 3.9, -0.1), (-0.1, 0.2, 4.2)]) * ang,
        (4, 4, 4),
        ['Ar', 'Ar'],
        [1.0, 1.0],
    )


def test_argon_crystal(argon_crystal):
    coords, lattice, k_grid, species, vol_ratios = argon_crystal()
    ene = MBDGeom(coords, lattice, k_grid, do_rpa=False).mbd_energy_species(species, vol_ratios, 0.83)
    assert ene == approx(-0.0021037562496878173, rel=1e-10)
    return ene

ene = test_argon_crystal(argon_crystal)

mbd_coords, mbd_lattice, k_grid, species, vol_ratios = argon_crystal()

print(mbd_coords)
print(mbd_lattice)

print(ene)
exit()
print(ang)
print(P_Ang)

print(species)
print(vol_ratios)

print(MBDGeom(atomic_positions, lattice, k_grid, do_rpa=False).mbd_energy_species(species, vol_ratios, 0.83))

exit()




