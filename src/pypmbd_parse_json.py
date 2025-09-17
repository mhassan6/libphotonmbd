import numpy as np
import json
from atomic_constants import P_Har, P_Ang
from photon_mode import photon_mode
from pypmbd_lib import xc_pmbd
from vdw.modelmbd import DEFAULT_RSSCS_BETA


def pypmbd_parse_json(filename, pt, gauge, noint = False):

   json_file = filename

   with open(json_file, "r") as f:
     data = json.load(f)

   code = data['code']
   xc_fun = data['xc']
   print(f'code: {code} functional: {xc_fun} mbd_beta: {DEFAULT_RSSCS_BETA[xc_fun]}')

   atomic_positions = P_Ang*np.array(data['atomic_positions']) #atomic_positions in a.u. while json in AA
   atomic_species = data['atomic_species']
   natm = len(atomic_positions)
   print('natm', natm)
   #atomic_species *= natm
   lattice = data['lattice'] #FIX to a.u.
   if lattice == 0:
      lattice = None
   else:
      lattice = P_Ang*np.array(lattice) #lattice in a.u. while json in AA
   kpoints = data['kpoints']
   if kpoints == 0:
      kpoints = None
   E0 = data['total_energy']
   xc_fun = data['xc']

   ratio = data['ratio'] # if using VASP
    
   print('atomic_positions: ', atomic_positions)
   print('atomic_species: ', atomic_species)
   print('lattice: ', lattice)
   print('kpoints: ', kpoints)
   print('Veff/Vfree ratio: ', ratio)
   print('npt modes:', pt.nmodes)
   print('pt freq:',pt.omega)
   print('pt lam:', pt.lam,'\n\n')
    
   output_array = xc_pmbd(atomic_positions, atomic_species, ratio, pt, beta = DEFAULT_RSSCS_BETA[xc_fun], gauge = gauge, \
      lattice=lattice, k_grid=kpoints, noint = noint)
   output_array.update({'E0': E0})
   output_array.update({'json': data})
   return output_array


#gauge = 'velocity'
#vdwene_len = xc_pmbd(atomic_positions, atomic_species, ratio, pt, beta = DEFAULT_RSSCS_BETA[xc_fun], gauge = gauge, lattice=lattice, k_grid=kpoints)
#print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)
