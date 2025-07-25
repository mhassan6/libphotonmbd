import numpy as np

class photon_mode:
  def __init__(self, nmodes, omega, lam, pol_vec, exc = 0., number = 0., qa = 0, pa = 0):
    self.nmodes = nmodes
    self.omega = omega
    self.lam = lam
    self.pol_vec = np.array(pol_vec)
    self.exc = exc
    self.number = number
    self.qa = qa #photon displacement
    self.pa = pa #photon momentum
