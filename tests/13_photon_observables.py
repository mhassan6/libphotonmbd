
import sys
sys.path.append("../src")
import numpy as np

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode

from pytest import approx

pt = photon_mode(1, 0.2, 0.1, [0,0,1])
filename = 'json/ar_dummy_natm_1.json'
gauge = 'velocity'
output_array = pypmbd_parse_json(filename, pt, gauge)

vdwene_len = output_array['ene']
E0 = output_array['E0']
pt_AT_energy = output_array['pt_AT_energy']
pt_rpa_orders = output_array['pt_rpa_orders']
pt_number = output_array['pt_number']
pt_adadaa = output_array['pt_adadaa']
pt_xx = output_array['pt_xx']
pt_pp = output_array['pt_pp']

print('\n')
print('with photons')
print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)

print('photonnumber', pt_number)
print('photon adadaa', pt_adadaa)

print('<x^2>', pt_xx)
print('<p^2>', pt_pp)

print('\n')
print('Loading data from mbdmatrix')
print('\n')

data = np.loadtxt('matrix/mbdmatrix.dat',delimiter=',')
#print(data)
data = np.delete(data, 0, axis=0)
data = np.delete(data, 0, axis=1)

data = np.delete(data, 0, axis=0)
data = np.delete(data, 0, axis=1)

data = np.delete(data, 2, axis=0)
data = np.delete(data, 2, axis=1)

data = np.delete(data, 2, axis=0)
data = np.delete(data, 2, axis=1)
#print(data)

from scipy import sparse as sp
import scipy.sparse.linalg as spla

def build_harmonic_oscillator_operators(N_pt):
  from scipy import sparse as sp
  import numpy as np
  import time

  unit = sp.identity(N_pt)
  create = np.zeros((N_pt,N_pt))
  destroy = np.zeros((N_pt,N_pt))

  
  for ii in range(1,N_pt):
    for jj in range(1,N_pt):
      if ii==jj:
        if ii <= N_pt-1:
          create[ii,jj-1]=np.sqrt(ii)

  for ii in range(0,N_pt):
    for jj in range(0,N_pt):
      if ii==jj:
        if jj <= N_pt-1:
          destroy[ii-1,jj]=np.sqrt(ii)

  return unit, sp.csr_matrix(create), sp.csr_matrix(destroy)

unit, create, destroy = build_harmonic_oscillator_operators(20)

create1 = sp.kron(unit, create)
create2 = sp.kron(create, unit)
destroy1 = sp.kron(unit, destroy)
destroy2 = sp.kron(destroy,unit)

oma = np.sqrt(data[2,2])
omp = pt.omega[0]

p1 = 1.0j*np.sqrt(oma/2)*(create1 - destroy1)
p2 = 1.0j*np.sqrt(omp/2)*(create2 - destroy2)

q1 = 1.0/np.sqrt(2*oma)*(create1 + destroy1)
q2 = 1.0/np.sqrt(2*omp)*(create2 + destroy2)

ham_int = 1/2*(data[0,0]*p1 @ p1 + data[1,1]*p2 @ p2 + data[0,1]*p1 @ p2 + data[1,0] * p2 @ p1) + 1/2*(data[2,2]*q1 @ q1 + data[3,3]*q2 @ q2)
ham_noint = 1/2*(p1 @ p1 + p2 @ p2) + 1/2*(data[2,2]*q1 @ q1 + data[3,3]*q2 @ q2)

eigs_noint, _ = np.linalg.eigh(ham_noint.todense())
eigs, evec = np.linalg.eigh(ham_int.todense())

print('energy',  pt.exc-vdwene_len, eigs[0] - eigs_noint[0])
print('photonnumber', pt_number, evec[:,0].T @ create2 @ destroy2 @ evec[:,0])
print('photon adadaa', pt_adadaa, evec[:,0].T @ create2 @ create2 @ destroy2 @ destroy2 @ evec[:,0])

print('<x^2>', pt_xx, evec[:,0].T @ q2 @ q2 @ evec[:,0])
print('<p^2>', pt_pp, evec[:,0].T @ p2 @ p2 @ evec[:,0])
print('<x*p>', evec[:,0].T @ p2 @ q2 @ evec[:,0])
print('<p*x>', evec[:,0].T @ q2 @ p2 @ evec[:,0])

print('Check for higher order correlations')
print('x**4', evec[:,0].T @ q2 @ q2 @ q2 @ q2 @ evec[:,0])
print('Gaussian: 3* <x**2>**2', 3*(evec[:,0].T @ q2 @ q2 @ evec[:,0])**2)

print('\nCheck for entropy \n')
rho = evec[:,0] @ evec[:,0].T

#print(rho.shape)
rho = np.array(rho)
drho = rho.reshape([20,20,20,20])
#print(drho.shape)
#drho =  np.einsum('ijkl->ik', drho)
drho =  np.einsum('jiki->jk', drho)

#print(drho.shape)

eigs_rho, _ = np.linalg.eigh(drho)

#print(eigs)

S1 = 0
for ii in range(len(eigs_rho)):
  if (eigs_rho[ii]) > 10**(-14):
    S1 += -eigs_rho[ii]*np.log(eigs_rho[ii])

nu = 2*np.sqrt(pt_xx*pt_pp)
S2 = (nu + 1) / 2 * np.log((nu + 1) / 2) -  (nu - 1) / 2 * np.log((nu - 1) / 2)


print('entropy', S1, S2)
print('')

assert pt.exc-vdwene_len == approx(eigs[0] - eigs_noint[0], abs=1e-10)
assert pt_number == approx(evec[:,0].T @ create2 @ destroy2 @ evec[:,0], abs=1e-10)
assert pt_adadaa == approx(evec[:,0].T @ create2 @ create2 @ destroy2 @ destroy2 @ evec[:,0], abs=1e-10)
assert pt_xx == approx(evec[:,0].T @ q2 @ q2 @ evec[:,0], abs=1e-10)
assert pt_pp == approx(evec[:,0].T @ p2 @ p2 @ evec[:,0], abs=1e-10)
assert S1 == approx(S2, abs = 1e-10)

log_val = (pt_xx * pt.omega) / (pt_pp / pt.omega)
r = 1/4 * np.log(log_val)

print('sqeezing: log_val', log_val, 'r', r)
print('')

print('13_photon_observables.py: Photon observables in velocity gauge tests passed!')

