
import sys
sys.path.append("../src")
import numpy as np

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode

from pytest import approx

filename = 'json/dummy_C_natm_8.json'
pt_2 = photon_mode(1, 0.2, 0.1/np.sqrt(2), [0,0,1])


gauge = 'length' #'velocity' #length'
output_array = pypmbd_parse_json(filename, pt_2, gauge,shift = 0)

vdwene_len2 = output_array['ene']
E0 = output_array['E0']
AT_energy = output_array['AT_energy']
rpa_orders = output_array['rpa_orders']
pt_AT_energy = output_array['pt_AT_energy']
pt_rpa_orders = output_array['pt_rpa_orders']
pt_number = output_array['pt_number']
pt_adadaa = output_array['pt_adadaa']

print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len2, 'vdw+pt', pt_2.exc)#, 'pt', pt.exc-vdwene_len)
print('photonnumber', pt_number)
print('photon adadaa', pt_adadaa)
print('')


pt_1 = photon_mode(1, 0.2, 0.1, [0,0,1])

filename = 'json/dummy_C_natm_4.json'
gauge = 'length' #'velocity' #'length'
output_array = pypmbd_parse_json(filename, pt_1, gauge, shift = 0)

vdwene_len1 = output_array['ene']
E0 = output_array['E0']
AT_energy = output_array['AT_energy']
rpa_orders = output_array['rpa_orders']
pt_AT_energy = output_array['pt_AT_energy']
pt_rpa_orders = output_array['pt_rpa_orders']
pt_number = output_array['pt_number']
pt_adadaa = output_array['pt_adadaa']

print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len1, 'vdw+pt', pt_1.exc)#, 'pt', pt.exc-vdwene_len)
print('photonnumber', pt_number)
print('photon adadaa', pt_adadaa)

assert vdwene_len1 == approx(vdwene_len2 / 2, abs=1e-10)

assert pt_1.exc == approx(pt_2.exc / 2, abs=1e-10)

print()
print('08_periodic_C_supercell_to_unitcell.py: scaling tests passed!')


