
import sys
sys.path.append("../src")
import numpy as np

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode
from atomic_constants import P_Har

pt = photon_mode(1, 2/P_Har, 0.05, [0,0,1])
filename = '../json/vasp_input_output_supercell6.json'
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

print(f'VASP:   Edisp (eV):   -6.18038, {-6.18038 / P_Har} (au)')

#atomic_positions = P_Ang*np.array(output_array['json']['atomic_positions']
#atomic_species = data['atomic_species'] #Fix to array
#natm = len(atomic_positions)
#atomic_species *= natm

#from pymbd.fortran import  MBDGeom
#geom = MBDGeom(atomic_positions, do_rpa=True, get_rpa_orders=True)
#ene, orders = geom.mbd_energy_species(atomic_species, ratio, DEFAULT_RSSCS_BETA[xc_fun])
#print('rsscs', ene, orders[1], orders[2], orders[3])
