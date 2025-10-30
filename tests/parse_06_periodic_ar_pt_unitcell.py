
import sys
sys.path.append("../src")
sys.path.append("../src/parse")
import numpy as np
import json
from atomic_constants import P_Ang


from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode


from pytest import approx
# omega=[0.1, 0.5,0.5]
# lam = [0.05, 0.05,0.05]
# pol = np.array([[0,0,1],[0,0,1],[0,0,1]]).T
# pt  = photon_mode(3, omega, lam, pol)


photonmodes = 'photonpilotexamples/3_modes'


pt = photon_mode.parse(photonmodes)



filename = 'json/ar_dummy_periodic_unitcell.json'
#filename = 'json/dummy_Ar_natm_2_dist_4_2k.json'
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


omega=[0.1, 0.5,0.5]
lam = [0.05, 0.05,0.05]
pol = np.array([[0,0,1],[0,0,1],[0,0,1]]).T
pt  = photon_mode(3, omega, lam, pol)




output_array1 = pypmbd_parse_json(filename, pt, gauge)


vdwene_len1 = output_array['ene']
E01 = output_array['E0']
AT_energy1 = output_array['AT_energy']
rpa_orders1 = output_array['rpa_orders']
pt_AT_energy1 = output_array['pt_AT_energy']
pt_rpa_orders1 = output_array['pt_rpa_orders']
pt_number1 = output_array['pt_number']
pt_adadaa1 = output_array['pt_adadaa']


print('Gauge', gauge, 'E0', E01, 'vdwene', vdwene_len1, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len1)


print('ene', vdwene_len, 'AT ene', AT_energy, 'ratio', vdwene_len1/AT_energy1)
print('log AT ene', AT_energy1, 'rpa_order', np.sum(rpa_orders1))


print('\n')
print('with photons')
print('ene', pt.exc, 'AT ene', pt_AT_energy1, 'ratio', pt.exc/pt_AT_energy1)
print('log AT ene', pt_AT_energy1, 'rpa_order', np.sum(pt_rpa_orders1))
print('photonnumber', pt_number1)
print('photon adadaa', pt_adadaa1)




print(vdwene_len1==vdwene_len)
print(E01==E0)
print(AT_energy1==AT_energy)
print(rpa_orders1==rpa_orders1)
print(pt_AT_energy1==pt_AT_energy)
print(pt_rpa_orders1==pt_rpa_orders)
print(pt_number1==pt_number)
print(pt_adadaa1==pt_adadaa)


print()


exit()
