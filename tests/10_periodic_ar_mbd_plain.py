
import sys
sys.path.append("../src")
import numpy as np
import json
from atomic_constants import P_Ang

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode

from pytest import approx

pt = photon_mode(1, 0.2, 0.1, [0,0,1])

filename = 'json/ar_dummy_periodic_unitcell.json'
#filename = 'json/dummy_Ar_natm_2_dist_4_2k.json'
#filename = 'json/dummy_Ar_natm_1_dist_4_4k.json'

gauge = 'length'
output_array = pypmbd_parse_json(filename, pt, gauge, shift = 0.5, runmode = 'plain')

vdwene_len = output_array['ene']
E0 = output_array['E0']
AT_energy = output_array['AT_energy']
rpa_orders = output_array['rpa_orders']
pt_AT_energy = output_array['pt_AT_energy']
pt_rpa_orders = output_array['pt_rpa_orders']
pt_number = output_array['pt_number']
pt_adadaa = output_array['pt_adadaa']

print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)

#Compare first to MBD

from pymbd.fortran import MBDGeom

with open(filename, "r") as f:
  data = json.load(f)
  atomic_positions = P_Ang*np.array(data['atomic_positions']) #atomic_positions in a.u. while json in AA
  lattice = data['lattice'] #FIX to a.u.
  lattice = P_Ang*np.array(lattice) #lattice in a.u. while json in AA
  kpoints = data['kpoints']

geom = MBDGeom(atomic_positions, lattice, kpoints)
mbd_ene = geom.mbd_energy_species(['Ar', 'Ar'], [1, 1], 0.83, variant = 'plain')

print(mbd_ene)
print(vdwene_len)
assert vdwene_len == approx(mbd_ene, rel=1e-10)

print()
print('10_periodic_ar_mbd_plain.py: MBD gauge tests passed!')

exit()



