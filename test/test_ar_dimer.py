
import sys
sys.path.append("../src")
import numpy as np

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode

pt = photon_mode(1, 0.2, 0.1, [0,0,1])
filename = '../json/ar_pyscf.json'
gauge = 'length'
output_array = pypmbd_parse_json(filename, pt, gauge)

vdwene_len = output_array['ene']
E0 = output_array['E0']
AT_energy = output_array['AT_energy']
rpa_orders = output_array['rpa_orders']
pt_AT_energy = output_array['pt_AT_energy']
pt_rpa_orders = output_array['pt_rpa_orders']

print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)

print('ene', vdwene_len, 'AT ene', AT_energy, 'ratio', vdwene_len/AT_energy)
print('log AT ene', AT_energy, 'rpa_order', np.sum(rpa_orders))

print('\n')
print('with photons')
print('ene', pt.exc, 'AT ene', pt_AT_energy, 'ratio', pt.exc/pt_AT_energy)
print('log AT ene', pt_AT_energy, 'rpa_order', np.sum(pt_rpa_orders))
