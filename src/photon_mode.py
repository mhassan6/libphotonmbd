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

  @classmethod
  def parse(self, filename):
    e = np.loadtxt(filename,max_rows=1, dtype=int)
    nmodes = int(e[0])
    cols = int(e[1])
    data = np.loadtxt(filename, skiprows=1, dtype=float)
    if cols != 5:
      print('Each mode requires frequency, wavelength, and X,Y,Z polarization')
    omega= [None]*nmodes
    lam= [None]*nmodes
    pol=np.empty((nmodes,3),int)
    for ii in range(0,nmodes):
      omega[ii]=data[ii,0]
      lam[ii]=data[ii,1]
      pol[ii,0]=data[ii,2]
      pol[ii,1]=data[ii,3]
      pol[ii,2]=data[ii,4]
    pol_vec=pol.T
    self.nmodes = nmodes
    self.omega = np.array(omega)
    self.lam = np.array(lam)
    self.pol_vec = np.array(pol)

  def __str__(self):
    return f"N Photonmodes: {self.nmodes}\n Omega: {self.omega} \n Lambda: {self.lam}"


