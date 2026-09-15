import json
import time
import numpy as np

from atomic_constants import P_Ang
from pypmbd_lib import xc_pmbd
from helper_functions import format_time, DEFAULT_RSSCS_BETA

#There are currently three rundmodes:
#rsscs, plain, plain2 (plain + dip,gg in dipole_matrix, see 10.1126/science.aae0509)

def pypmbd_parse_json(filename, pt, gauge, noint = False, rpa = True, runmode = 'rsscs', shift = 0., custom_kpts = None):

   start_time = time.time()

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
   elif custom_kpts is not None:
      kpoints = custom_kpts
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
   print('pt lam:', pt.lam)
   print('shift:', shift, '\n\n')
    
   output_array = xc_pmbd(atomic_positions, atomic_species, ratio, pt, beta = DEFAULT_RSSCS_BETA[xc_fun], gauge = gauge, \
      lattice=lattice, k_grid=kpoints, noint = noint, rpa = rpa, runmode = runmode, shift = shift)
   output_array.update({'E0': E0})
   output_array.update({'json': data})

   print('')
   elapsed_seconds = time.time() - start_time
   print('PMBD Runtime: --- ' + format_time(elapsed_seconds) + ' hours:minutes:seconds ---')
   print("PMBD Runtime: --- %s seconds ---" %elapsed_seconds)
   print('')

   return output_array


#gauge = 'velocity'
#vdwene_len = xc_pmbd(atomic_positions, atomic_species, ratio, pt, beta = DEFAULT_RSSCS_BETA[xc_fun], gauge = gauge, lattice=lattice, k_grid=kpoints)
#print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)
