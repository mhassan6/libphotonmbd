import sys
sys.path.append("../src")
sys.path.append("../src/parse")
import numpy as np
import json
from atomic_constants import P_Ang

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode

#
#0.134 0.04435 1 0 0
#0.52 0.02123 0 1 0
#0.81 0.056547 0 0 1

omega=[0.134, 0.52,0.81]
lam = [0.05, 0.05,0.05]
pol = np.array([[0,0,1],[0,0,1],[0,0,1]]).T
pt1  = photon_mode(3, omega, lam, pol)

photonmodes = 'photonmodes/three_modes.dat'

pt2 = photon_mode(photonmodes)

assert pt1.nmodes==pt2.nmodes
assert pt1.omega==pt2.omega

print('01_test_finite_ar_dimer.py: Velocity gauge tests passed!')


#print(pt1.lam==pt2.lam)
#print(pt1.pol_vec==pt2.pol_vec)
#print(pt1.nmodes==pt2.nmodes)
#print(pt1.exc==pt2.exc)
#print(pt1.number==pt2.number)
#print(pt1.qa==pt2.qa)
#print(pt1.pa==pt2.pa)
