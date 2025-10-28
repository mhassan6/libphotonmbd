import sys
sys.path.append("../src")
import numpy as np

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode

from pytest import approx

gauge = 'length' #'velocity' #length'


for ii in range(0,3):
  if ii == 0:
    pol_vec = [0,0,1]
  elif ii == 1:
    pol_vec = [1,0,0]
  elif ii == 2:
    pol_vec = [0,1,0]

  filename = 'json/dummy_C_natm_8.json'
  pt_2 = photon_mode(1, 0.2, 0.1/np.sqrt(2), pol_vec)
  output_array = pypmbd_parse_json(filename, pt_2, gauge,shift = 0)

  vdwene_len2 = output_array['ene']
  E0 = output_array['E0']

  print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len2, 'vdw+pt', pt_2.exc)#, 'pt', pt.exc-vdwene_len)
  print('')

  pt_1 = photon_mode(1, 0.2, 0.1, pol_vec)

  filename = 'json/dummy_C_natm_4.json'
  output_array = pypmbd_parse_json(filename, pt_1, gauge, shift = 0)

  vdwene_len1 = output_array['ene']
  E0 = output_array['E0']

  print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len1, 'vdw+pt', pt_1.exc)#, 'pt', pt.exc-vdwene_len)
  print('')

  assert vdwene_len1 == approx(vdwene_len2 / 2, abs=1e-10)
  assert pt_1.exc == approx(pt_2.exc / 2, abs=1e-10)

  print()
  print(f'08_periodic_C_supercell_to_unitcell.py: scaling test {ii}  passed!')


