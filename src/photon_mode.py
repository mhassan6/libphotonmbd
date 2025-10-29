import numpy as np

class photon_mode:
  def __init__(self, nmodes, omega, lam, pol_vec, exc = 0., number = 0., qa = 0, pa = 0):
    self.nmodes = nmodes
    if nmodes == 1:
      self.omega = np.array([omega])
      self.lam = np.array([lam])
      #self.pol_vec = pol_vec
      self.pol_vec = np.array([pol_vec]).T
    else:
      self.omega = omega
      self.lam = lam
      self.pol_vec = np.array(pol_vec)

    self.exc = exc
    self.number = number
    self.qa = qa #photon displacement
    self.pa = pa #photon momentum
  def parse(filename):
    d = np.loadtxt(filename, skiprows=1, dtype=float)
    e=np.loadtxt(filename,max_rows=1, dtype=int)
    N=int(e[0])
    j=int(e[1])
    if j!=5:
        print('Each mode requires frequency, wavelength, and X,Y,Z polarization')

    omega= [None]*N
    lam= [None]*N
    pol=np.empty((N,3),int)
    for i in range(0,N):
        omega[i]=d[i,0]
        lam[i]=d[i,1]
        pol[i,0]=d[i,2]
        pol[i,1]=d[i,3]
        pol[i,2]=d[i,4]
    pol_vec=pol.T
    return photon_mode(N, omega, lam, pol_vec)
