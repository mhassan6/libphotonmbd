#compare to one test from here:
#https://github.com/libmbd/libmbd/blob/master/tests/test_fortran.py

#def test_argon_dimer_rsscs_rpa():
#    geom = MBDGeom([(0, 0, 0), (0, 0, 4 * ang)], do_rpa=True, get_rpa_orders=True)
#    ene, orders = geom.mbd_energy_species(['Ar', 'Ar'], [1, 1], 0.83)
#    assert ene == approx(-0.0002462647623815428, rel=1e-10)
#    assert orders[1] == approx(-0.0002461558113413099, rel=1e-10)
#    assert orders[2] == approx(0)
#    assert orders[3] == approx(-1.0885208380438466e-07, rel=1e-10)

import sys
sys.path.append("../src")
import numpy as np
import json
from atomic_constants import P_Ang

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode

from pytest import approx

pt = photon_mode(1, 0.2, 0.1, [0,0,1])
filename = 'json/ar_dummy.json'
gauge = 'length'
output_array = pypmbd_parse_json(filename, pt, gauge, rpa = True, runmode = 'plain')

vdwene_len = output_array['ene']
E0 = output_array['E0']
AT_energy = output_array['AT_energy']
pt_AT_energy = output_array['pt_AT_energy']
rpa_orders = output_array['rpa_orders']

print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)

assert vdwene_len == approx(AT_energy, rel=1e-10) #should be rel, but might be related to the ang factor
assert pt.exc == approx(pt_AT_energy, rel=1e-09) #should be rel, but might be related to the ang factor

from pymbd.fortran import MBDGeom

with open(filename, "r") as f:
  data = json.load(f)
  atomic_positions = P_Ang*np.array(data['atomic_positions']) #atomic_positions in a.u. while json in AA
  mbd_ene, mbd_orders = MBDGeom(atomic_positions, do_rpa=True, get_rpa_orders=True).mbd_energy_species(['Ar', 'Ar'], [1, 1], 0.83, variant='plain')


assert vdwene_len == approx(mbd_ene, rel=1e-10) #should be rel, but might be related to the ang factor
print('rpa_orders', rpa_orders)
assert rpa_orders[2]  == approx(mbd_orders[1], rel=1e-10)
assert rpa_orders[3]  == approx(mbd_orders[2], rel=1e-10)
assert rpa_orders[4]  == approx(mbd_orders[3], rel=1e-10)


#Also test plain2
output_array = pypmbd_parse_json(filename, pt, gauge, rpa = True, runmode = 'plain2')

vdwene_len = output_array['ene']
E0 = output_array['E0']
AT_energy = output_array['AT_energy']
pt_AT_energy = output_array['pt_AT_energy']

print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)

assert vdwene_len == approx(AT_energy, rel=1e-10) #should be rel, but might be related to the ang factor
assert pt.exc == approx(pt_AT_energy, rel=1e-09) #should be rel, but might be related to the ang factor

print()
print('12_finite_ar_dimer_mbd_plain_rpa.py: MBD tests passed!')

