import numpy as np
import json
from atomic_constants import P_Har, P_Ang
from photon_mode import photon_mode
from pypmbd_lib import xc_pmbd
from vdw.modelmbd import DEFAULT_RSSCS_BETA

json_file = '../json/vasp_input_output_supercell6.json'
#json_file = '../json/vasp_input_output_kpoint12.json'

with open(json_file, "r") as f:
  data = json.load(f)

xc_fun = 'pbe'
print(f'functional: {xc_fun} mbd_beta: {DEFAULT_RSSCS_BETA[xc_fun]}')

atomic_positions = P_Ang*np.array(data['atomic_positions']) #Fix to a.u.
atomic_species = data['atomic_species'] #Fix to array
natm = len(atomic_positions)
atomic_species *= natm
lattice = P_Ang*np.array(data['lattice']) #FIX to a.u.
kpoints = data['kpoints']
E0 = data['total_energy']/P_Har # if using VASP

ratio = data['ratio'] # if using VASP
    
print('atomic_positions: ', atomic_positions)
print('atomic_species: ', atomic_species)
print('lattice: ', lattice)
print('kpoints: ', kpoints)
print('Veff/Vfree ratio: ', ratio, '\n\n')
    
pt = photon_mode(1, 2/P_Har, 0.05, [0,0,1])

gauge = 'length'
vdwene_len = xc_pmbd(atomic_positions, atomic_species, ratio, pt, beta = DEFAULT_RSSCS_BETA[xc_fun], gauge = gauge, lattice=lattice, k_grid=kpoints)
print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)

#gauge = 'velocity'
#vdwene_len = xc_pmbd(atomic_positions, atomic_species, ratio, pt, beta = DEFAULT_RSSCS_BETA[xc_fun], gauge = gauge, lattice=lattice, k_grid=kpoints)
#print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)
