import sys
sys.path.append("../src")
sys.path.append("../src/parse")
import numpy as np
import json
from atomic_constants import P_Ang

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode



omega=[0.1, 0.5,0.5]
lam = [0.05, 0.05,0.05]
pol = np.array([[0,0,1],[0,0,1],[0,0,1]]).T
pt1  = photon_mode(3, omega, lam, pol)

photonmodes = 'photonpilotexamples/3_modes'

pt2 = photon_mode.parse(photonmodes)

print(pt1.nmodes==pt2.nmodes)
print(pt1.omega==pt2.omega)
print(pt1.lam==pt2.lam)
print(pt1.pol_vec==pt2.pol_vec)
print(pt1.nmodes==pt2.nmodes)
print(pt1.exc==pt2.exc)
print(pt1.number==pt2.number)
print(pt1.qa==pt2.qa)
print(pt1.pa==pt2.pa)
