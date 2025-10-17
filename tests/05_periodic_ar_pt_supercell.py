
import sys
sys.path.append("../src")
import numpy as np
import json
from atomic_constants import P_Ang

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode

from pytest import approx

pt = photon_mode(1, 0.2, 0.1, [0,0,1])

filename = 'json/ar_dummy_periodic_supercell.json'

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

print('')
print('')
pt_vel = photon_mode(1, 0.2, 0.1, [0,0,1])
gauge = 'velocity'
output_array = pypmbd_parse_json(filename, pt_vel, gauge)

vdwene_vel = output_array['ene']
E0 = output_array['E0']
AT_energy_vel = output_array['AT_energy']
rpa_orders = output_array['rpa_orders']
pt_AT_energy_vel = output_array['pt_AT_energy']
pt_rpa_orders = output_array['pt_rpa_orders']
pt_number = output_array['pt_number']
pt_adadaa = output_array['pt_adadaa']

print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)

print('ene', vdwene_vel, 'AT ene', AT_energy, 'ratio', vdwene_len/AT_energy)
print('log AT ene', AT_energy, 'rpa_order', np.sum(rpa_orders))

print('\n')
print('with photons')
print('ene', pt.exc, 'AT ene', pt_AT_energy, 'ratio', pt.exc/pt_AT_energy)
print('log AT ene', pt_AT_energy, 'rpa_order', np.sum(pt_rpa_orders))
print('photonnumber', pt_number)
print('photon adadaa', pt_adadaa)

#comparing length gauge vs velocity gauge values against each other
assert vdwene_len == approx(vdwene_vel, abs=1e-10)
assert pt.exc == approx(pt_vel.exc, abs=1e-10)

assert AT_energy == approx(AT_energy_vel, abs=1e-10)
assert pt_AT_energy == approx(pt_AT_energy_vel, abs=1e-10)

print()
print('05_periodic_ar_pt_supercell: Supercell@Gamma length gauge and Velocity gauge tests passed!')




