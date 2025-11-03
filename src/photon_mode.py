import numpy as np

class photon_mode:
  def __init__(self, *args, **kwargs): #nmodes, omega, lam, pol_vec):
    if len(args) == 1 and isinstance(args[0], str):
      print('Filename:', args[0])
      self.parse(args[0])
    else:
      self.nmodes = args[0]
      if self.nmodes == 1:
        self.omega = np.array([args[1]])
        self.lam = np.array([args[2]])
        #self.pol_vec = pol_vec
        self.pol_vec = np.array([args[3]]).T
      else:
        self.omega = args[1]
        self.lam = args[2]
        self.pol_vec = np.array(args[3])

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
