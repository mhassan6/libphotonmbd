import numpy as np
from pyscf import dft, gto, scf
from atomic_constants import P_Har, P_a_B
from vdw.modelmbd import DEFAULT_RSSCS_BETA
from hirshfeld import HirshfeldAnalysis

element = 'Xe'
atom_str = f'{element} 0 0 0; ' 
dist_ang = 4

for ii in range(1, 201): #(1, 20):
  dist = ii*dist_ang
  print('dist: ', dist)
  atom_str += f'{element} 0 0 {dist}; ' 
  print(atom_str)
  mol = gto.M(atom=atom_str)
  mol.basis = '3-21g' #ccpvdz' #augccpvdz'
  mol.build()
  print(mol._atom)
  mf = dft.RKS(mol) 

  xc_fun = 'pbe0'
  print(f'functional: {xc_fun} mbd_beta: {DEFAULT_RSSCS_BETA[xc_fun]}')

  #mf.xc = xc_fun

  E0 = 0 #mf.kernel()
  #hirshfeld_result = HirshfeldAnalysis(mf).run().result
  #Vfree, Veff = [hirshfeld_result[token] for token in ("V_free", "V_eff")]
  ratio = [1]*len(mf.mol.elements) #Veff/Vfree
  print('Veff/Vfree ratio', ratio)

  import json
  results = {
    "code": "dummy",
    "xc": xc_fun,
    "atomic_species": mf.mol.elements,
    "atomic_positions": mf.mol.atom_coords().tolist(),
    "lattice": 0,
    "kpoints": 0,
    "total_energy": E0,
    "ratio": ratio#.tolist(),
    }

  with open(f"json/dummy_{element}_natm_{ii+1}_dist_{dist_ang}.json", "w") as f:
        json.dump(results, f, indent=2)

  
