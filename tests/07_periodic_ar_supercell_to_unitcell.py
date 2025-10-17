
import sys
sys.path.append("../src")
import numpy as np

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode

from pytest import approx

filename = 'json/dummy_Ar_natm_4_dist_4_1k.json'
pt_4 = photon_mode(1, 0.2, 0.1/np.sqrt(4), [0,0,1])


gauge = 'length' #'velocity' #length'
output_array = pypmbd_parse_json(filename, pt_4, gauge, shift = 0)

vdwene_len4 = output_array['ene']
E0 = output_array['E0']
AT_energy = output_array['AT_energy']
rpa_orders = output_array['rpa_orders']
pt_AT_energy = output_array['pt_AT_energy']
pt_rpa_orders = output_array['pt_rpa_orders']
pt_number = output_array['pt_number']
pt_adadaa = output_array['pt_adadaa']

print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len4, 'vdw+pt', pt_4.exc)#, 'pt', pt.exc-vdwene_len)
print('photonnumber', pt_number)
print('photon adadaa', pt_adadaa)
print('')


filename = 'json/dummy_Ar_natm_2_dist_4_2k.json'
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

filename = 'json/dummy_Ar_natm_1_dist_4_4k.json'
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
assert vdwene_len1 == approx(vdwene_len4 / 4, abs=1e-10)

assert pt_1.exc == approx(pt_2.exc / 2, abs=1e-10)
assert pt_1.exc == approx(pt_4.exc / 4, abs=1e-10)

print()
print('07_periodic_ar_supercell_to_unitcell.py: scaling tests passed!')


