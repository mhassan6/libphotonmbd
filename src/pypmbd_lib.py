import numpy as np

from pymbd.pymbd import freq_grid, dipole_matrix, get_kpts, screening
from pymbd.fortran import MBDGeom

from pymbd import from_volumes#, mbd_energy
from pymbd_mod import mbd_energy_get_matrix2 #mbd_energy_get_matrix, screening_all_freqs, screening_all_freqs_orig

### finite difference:
def get_photon_number_fd(coords, elements, ratio, pt, beta):
  from photon_mode import photon_mode

  delta = 1e-3
  om = pt.omega
  lam = pt.lam

  pt_op = photon_mode(1, om + delta, lam, pt.pol_vec)
  vdwene = xc_mbd(coords, elements, ratio, pt_op, beta)
  pt_om = photon_mode(1, om - delta, lam, pt.pol_vec)
  vdwene = xc_mbd(coords, elements, ratio, pt_om, beta)

  pt_lp = photon_mode(1, om, lam + delta/10, pt.pol_vec)
  vdwene = xc_mbd(coords, elements, ratio, pt_lp, beta)

  pt_lm = photon_mode(1, om, lam - delta/10, pt.pol_vec)
  vdwene = xc_mbd(coords, elements, ratio, pt_lm, beta)

  ptnumber_mom = 0.5*(pt_op.exc - pt_om.exc)/delta + 0.25*lam/om*(pt_lp.exc - pt_lm.exc)/delta*10
  print('ptnumber in velocity gauge via finite difference', ptnumber_mom)
  return

### xc functions:

def get_pt_mbdmatrix(natm, pt, matrix, alpha, omega, gauge = 'length', noint = False, k_point_idx = 0):

  if noint: 
    matrix = np.diag(np.diag(matrix))

  if gauge == 'length':
      mbdmatrix = np.zeros((3*natm + pt.nmodes, 3*natm + pt.nmodes), dtype=complex)
      for ii in range(0, natm):
        for xx in range(0, 3):
          for jj in range(0, natm):
            for yy in range(0, 3):
              mbdmatrix[3*ii+xx,3*jj+yy] = matrix[3*ii+xx,3*jj+yy]
              #R2 term (add only for first kpoint (should be gamma?)
              if k_point_idx == 0:
                for imodes in range(0, pt.nmodes):
                  mbdmatrix[3*ii+xx,3*jj+yy] += pt.lam[imodes]**2*pt.pol_vec[xx,imodes]*pt.pol_vec[yy,imodes] \
                    *np.sqrt(alpha[ii])*np.sqrt(alpha[jj])*omega[ii]*omega[jj]

      for imodes in range(0, pt.nmodes):
        mbdmatrix[3*natm+imodes, 3*natm+imodes] = pt.omega[imodes]**2
      if k_point_idx == 0:
        for ii in range(0, natm):
          for xx in range(0, 3):
            for imodes in range(0, pt.nmodes):
              offdiagonal=-pt.lam[imodes]*pt.pol_vec[xx,imodes]*np.sqrt(alpha[ii])*omega[ii]*pt.omega[imodes]
              mbdmatrix[3*natm+imodes, 3*ii+xx] += offdiagonal
              mbdmatrix[3*ii+xx, 3*natm+imodes] += offdiagonal
  
  elif gauge == 'velocity':
    mbdmatrix = np.zeros((2*(3*natm + pt.nmodes), 2*(3*natm + pt.nmodes)), dtype = complex)
    mbdmatrix[(3*natm + pt.nmodes):(6*natm + pt.nmodes), (3*natm + pt.nmodes):(6*natm + pt.nmodes)] = matrix
    mbdmatrix[:(3*natm + pt.nmodes), :(3*natm + pt.nmodes)] = np.diag(np.ones((3*natm + pt.nmodes)))
    for imodes in range(0, pt.nmodes):
      mbdmatrix[6*natm + pt.nmodes + imodes, 6*natm + pt.nmodes + imodes] += pt.omega[imodes]**2

    if k_point_idx == 0:
      for ii in range(0, natm):
        for xx in range(0, 3):
          for imodes in range(0, pt.nmodes):
            offdiagonal = -pt.lam[imodes]*pt.pol_vec[xx,imodes]*np.sqrt(alpha[ii])*omega[ii]/pt.omega[imodes]
#           mbdmatrix[3*ii + xx + imodes, 3*natm + imodes] += offdiagonal
#           mbdmatrix[3*natm + imodes, 3*ii + xx + imodes] += offdiagonal
            mbdmatrix[3*ii + xx, 3*natm + imodes] += offdiagonal
            mbdmatrix[3*natm + imodes, 3*ii + xx] += offdiagonal

      for imodes in range(0, pt.nmodes):
        for ii in range(0, natm):
          for xx in range(0, 3):
            mbdmatrix[3*natm + imodes, 3*natm + imodes] += (pt.lam[imodes]*pt.pol_vec[xx,imodes])**2*alpha[ii]*omega[ii]**2/pt.omega[imodes]**2

    #for imodes in range(0, pt.nmodes):
    #  for jmodes in range(0,pt.nmodes):
    #    element = 0
    #    for ii in range(0, natm):
    #      for xx in range(0, 3):
    #        for jj in range(0, natm):
    #          for yy in range(0,3):  
    #            element += (pt.lam[imodes]*pt.pol_vec[xx,imodes])*(pt.lam[jmodes]*pt.pol_vec[yy,jmodes])*np.sqrt(alpha[ii])*np.sqrt(alpha[jj])*omega[ii]*omega[jj]/pt.omega[imodes]/pt.omega[jmodes]
    #    mbdmatrix[3*natm + imodes, 3*natm + jmodes] += element

  #np.savetxt(f"mbdmatrix.dat", np.real(mbdmatrix), delimiter=",", fmt='%.10f')

  return mbdmatrix

def get_pqtoada(nmodes):
  pqtoada = 1.0j*np.zeros((2*nmodes,2*nmodes)) #(a1, a2, ad1, ad2)
  for ii in range(0, nmodes):
    pqtoada[ii,ii] = 1./np.sqrt(2.) 
    pqtoada[ii,nmodes + ii] = 1./np.sqrt(2.)
    pqtoada[nmodes+ii,ii] = 1.0j*np.sqrt(1/2)
    pqtoada[nmodes+ii,nmodes+ii] = -1.0j*np.sqrt(1/2)
  return pqtoada

def diagonalize_pmbdmatrix(matrix, om):
  #J.H.P Colpa, Physica 93A 327-353 (1978)

  nmodes = int(len(matrix) / 2)
  #print('nmodes', nmodes)

  A = np.real(matrix[:nmodes,:nmodes])
  B = np.real(matrix[nmodes:,:nmodes])
  D = np.real(matrix)

  L = np.linalg.cholesky(A-B)
  #alternatively below via diagonalization
  #Lval, L = np.linalg.eig(A-B)
  #L = np.sqrt(Lval)*L

  # A - B = K K^dagger
  eval3, evec3 = np.linalg.eig(L.T.conj() @ (A+B) @ L)
  phi1 = 0j*evec3
  psi1 = 0j*evec3
  uvec = 0j*evec3
  vvec = 0j*evec3
  evec4 = 0j*evec3
  X_bo = 0j*evec3
  Y_bo = 0j*evec3

  for ii in range(0, nmodes):
    evec4[:,ii] = evec3[:,ii] / (eval3[ii])**(1/4)
    phi1[:,ii] = L @ evec4[:,ii]
    psi1[:,ii] = (A + B) @ phi1[:,ii] / np.sqrt(eval3[ii])
    uvec[:,ii] = 1/2*(phi1[:,ii] + psi1[:,ii])
    vvec[:,ii] = 1/2*(phi1[:,ii] - psi1[:,ii])


  Tm_new = 1.0j*np.zeros((2*nmodes,2*nmodes))
  Tm_new[:nmodes,:nmodes] = uvec
  Tm_new[nmodes:,:nmodes] = vvec
  Tm_new[:nmodes,nmodes:] = vvec
  Tm_new[nmodes:,nmodes:] = uvec

  evals = np.real(np.diag(Tm_new.T @ D @ Tm_new))

  X_bo = 1/2*(phi1.T @ np.diag(om**(-0.5)) + psi1.T @ np.diag(om**(0.5)))
  Y_bo = 1/2*(phi1.T @ np.diag(om**(-0.5)) - psi1.T @ np.diag(om**(0.5)))

  return evals, np.real(phi1.T), np.real(psi1.T)

def diagonalize_pmbdmatrix_gamma(matrix, om):

  nmodes = int(len(matrix) / 2)
  gamma = np.zeros((2*nmodes, 2*nmodes))
  gamma[:nmodes, :nmodes] = np.diag(np.ones(nmodes))
  gamma[nmodes:, nmodes:] = -np.diag(np.ones(nmodes))

  eigvalada, eigvecada = np.linalg.eig(np.real(gamma @ matrix))
  evals = abs(eigvalada)
  return evals, 0*gamma[:,nmodes].T, 0*gamma[:,nmodes].T

def get_omegas(omega_rsscs, ptomega, natm, ptnmodes):
    #create omegas array for photon number
    omegas = np.zeros((3*natm + ptnmodes,))
    count = 0
    for ii in range(0, natm):
      for xx in range(0,3):
          omegas[count] = omega_rsscs[ii]
          count += 1
    omegas[3*natm:] = ptomega
    return omegas


def calc_ptobservables(omegas, phi1, psi1, natm):
    #this assumes the photon number is in the last column/row
    #ombar = np.sqrt(eigvalq) #*2

    #Y1 = np.diag(ombar**(0.5)) @ eigvecq.T @ np.diag(omegas**(-0.5))
    #Y2 = np.diag(ombar**(-0.5)) @ eigvecq.T @ np.diag(omegas**(0.5))

    Y1 = phi1 @ np.diag(omegas**(-0.5))
    Y2 = psi1 @ np.diag(omegas**0.5)

    X = (Y1 + Y2) / 2
    Y = (Y1 - Y2) / 2

    ptnumber_matrix = Y.T @ Y
    ptnumber = np.sum(np.diag(ptnumber_matrix[3*natm:, 3*natm:]))
    adadaambd=(Y.T@ X)[-1,-1]*(X.T @ Y)[-1,-1]+2*(Y.T @ Y)[-1,-1]**2 #+(Y.T @ Y)[-1,-1]* (Y @ Y)[-1,-1]

    YYt = Y.T @ Y
    XtY = X.T @ Y
    xx = 1/(omegas[-1])*(YYt[-1,-1]-XtY[-1,-1]+1/2) #careful assumes only one photon mode
    pp = (omegas[-1])*(YYt[-1,-1]+XtY[-1,-1]+1/2)   #careful assumes only one photon mode

    return ptnumber, adadaambd, X, Y, xx, pp #[-1, -1]

def get_AT_matrix(coords, alpha_0, C6, R_vdw, beta, pt = None, lattice = None, k_grid = None, noint = False, runmode = 'rsscs', shift = 0.5):
    #see https://github.com/libmbd/libmbd/blob/master/src/mbd_rpa.F90
    nfreq = 15
    natm = alpha_0.shape[0]
    rpa_order_max = 10
    rpa_orders = np.zeros(rpa_order_max)

    freq, freq_w = freq_grid(nfreq)

    if runmode == 'rsscs':
      alpha_0_rsscs, C6_rsscs, R_vdw_rsscs = screening(
        coords, alpha_0, C6, R_vdw, beta, lattice=lattice, nfreq=15)
      omega_rsscs = 4 / 3 * C6_rsscs / alpha_0_rsscs**2
    elif runmode == 'plain' or runmode == 'plain2':
      alpha_0_rsscs = alpha_0
      omega_rsscs = 4 / 3 * C6 / alpha_0**2
      R_vdw_rsscs = R_vdw
    else:
      print('This runmode is not supported!')
      exit()

#    alpha_0_rsscs, C6_rsscs, R_vdw_rsscs = screening(
#        coords, alpha_0, C6, R_vdw, beta, lattice = lattice, nfreq=nfreq
#    )
#    omega = 4 / 3 * C6_rsscs / alpha_0_rsscs**2
    
    alpha_dyn = [alpha_0_rsscs / (1 + (u / omega_rsscs) ** 2) for u in freq]
    
    if lattice is None:
        k_points = [None]
        Mgeom = MBDGeom(coords)
    else:
        assert k_grid is not None
        k_points = get_kpts(lattice, k_grid, shift = shift)
        Mgeom = MBDGeom(coords, lattice, k_points)

    if runmode == 'plain2':
      sigma_array = []
      for a in alpha_0_rsscs:
        sigma = (np.sqrt(2 / np.pi) * a / 3) ** (1 / 3)
        sigma_array.append(sigma)

    ene = 0.0
    for idx, k_point in enumerate(k_points):
#    for k_point in k_points:
      if not noint:
        if runmode == 'plain2':
          dipole_matrix_mbd = Mgeom.dipole_matrix('dip,gg', sigma=sigma_array, k_point=k_point)
        else:
         dipole_matrix_mbd = Mgeom.dipole_matrix('fermi,dip', beta=beta, k_point=k_point, R_vdw= R_vdw_rsscs)

#        dipole_matrix_mbd = dipole_matrix(
#                  coords,
#                  'fermi,dip',
#                  R_vdw=R_vdw_rsscs,
#                  beta=beta,
#                  lattice=lattice,
#                  k_point=k_point)
      else:
        dipole_matrix_mbd = np.zeros((3*natm, 3*natm))

      energy = 0.0

      for ifreq in range(1, nfreq + 1):
        #this is what is in e.g. Bucko et al 2016 J.Phys.:Condens.Matter 28 045201
        #AT = np.diag(np.repeat(-alpha_0_rsscs[ifreq], 3)) @ dipole_matrix

        #this seems to be what is actually implemented in mbd_rpa (see above)
        
        AT = np.zeros((3*natm, 3*natm), dtype = complex)
        for ii in range(0, natm):
          for xx in range(0, 3):
            for jj in range(0, natm):
              for yy in range(0, 3):
                res = 0  
                if pt is not None and idx == 0:
                  for imodes in range(0, pt.nmodes):
                    res += pt.lam[imodes]**2*pt.pol_vec[xx,imodes]*pt.pol_vec[yy,imodes]*freq[ifreq]**2/(freq[ifreq]**2 + pt.omega[imodes]**2)
                AT[3*ii+xx,3*jj+yy] = np.sqrt(alpha_dyn[ifreq][ii]*alpha_dyn[ifreq][jj]) \
                    *(dipole_matrix_mbd[3*ii+xx,3*jj+yy] + res)

        eigs = np.linalg.eigvalsh(AT)
        log_eigs = np.log(1.0 + eigs)
        energy += (1/2/np.pi)*np.sum(log_eigs)*freq_w[ifreq]

        for irpa in range(1, rpa_order_max):
          rpa_orders[irpa] += -1/2/np.pi*(-1)**irpa*np.sum(eigs**irpa / irpa)*freq_w[ifreq]  #fix /len(k_points)

      ene += energy / len(k_points)

    rpa_orders = rpa_orders / len(k_points)
      
    return ene, rpa_orders

def xc_pmbd(coords, elements, ratio, pt = None, beta = 0.84, gauge = 'length', lattice=None, k_grid=None, noint = False, rpa = True, runmode = 'rsscs', shift=0.5): #gauge = length, velocity

    natm = len(elements)
    alpha_0, C6, R_vdw = from_volumes(elements, ratio)

    if lattice is None:
        k_points = [None]
    else:
        assert k_grid is not None
        k_points = get_kpts(lattice, k_grid, shift = shift)
        print('k_points: ', k_points)
        print('len k_points: ', len(k_points))
 
    ene, alpha_0_rsscs, omega_rsscs, matrix, eigvecs_list, eigs_list = mbd_energy_get_matrix2(coords, alpha_0, C6, R_vdw, \
                                                 beta=beta, lattice=lattice, k_grid=k_grid, noint = noint, \
                                                 runmode = runmode, shift = shift)

    if runmode == 'plain' or runmode == 'plain2':
      omega_rsscs = 4 / 3 * C6 / alpha_0**2
      alpha_0_rsscs = alpha_0

    if rpa is True:
      AT_energy, rpa_orders = get_AT_matrix(coords, alpha_0, C6, R_vdw, beta=beta, lattice=lattice, k_grid=k_grid, noint = noint, runmode = runmode, shift = shift)
    else:
      AT_energy = 0
      rpa_orders = 0

    omegas = get_omegas(omega_rsscs, pt.omega, natm, pt.nmodes)
    pt.exc = 0.0
    psi1 = np.zeros((3*natm + pt.nmodes, 3*natm + pt.nmodes), dtype = complex)
    phi1 = np.zeros((3*natm + pt.nmodes, 3*natm + pt.nmodes), dtype = complex)

    if gauge == 'length':
      for idx in range(0, len(k_points)):
        if idx == 0:
          mbdmatrix = get_pt_mbdmatrix(natm, pt, matrix[idx], alpha_0_rsscs, omega_rsscs, gauge, noint, idx)
          eigs, eigvecs = np.linalg.eigh(mbdmatrix)
        else:
          eigvecs = np.zeros((3*natm + pt.nmodes, 3*natm + pt.nmodes), dtype = complex)
          eigvecs[:3*natm, :3*natm] = eigvecs_list[idx]
          eigvecs[3*natm:, 3*natm:] = np.diag(np.ones(pt.nmodes))

          eigs = np.zeros((3*natm + pt.nmodes))
          eigs[:3*natm] = eigs_list[idx]
          eigs[3*natm:] = [x**2 for x in pt.omega]
 
        pt.exc += np.sum(np.sqrt(eigs)) / 2

        #print(idx, 'eigs', np.sqrt(eigs))
        
        ombar = np.sqrt(eigs) #*2
        phi1 += np.diag(ombar**(0.5)) @ eigvecs.T
        psi1 += np.diag(ombar**(-0.5)) @ eigvecs.T

    elif gauge == 'velocity':
      if pt.nmodes > 1:
        print('more than 1 photon mode is not implemented yet for velocity gauge')
        exit()
      pqtoada = get_pqtoada(3*natm + pt.nmodes)
      for idx in range(0, len(k_points)): 
        mbdmatrix = get_pt_mbdmatrix(natm, pt, matrix[idx], alpha_0_rsscs, omega_rsscs, gauge, noint, idx)

        eigs, phi1_t, psi1_t = diagonalize_pmbdmatrix(np.conj(pqtoada).T @ mbdmatrix @ pqtoada, omegas)
        #eigs, phi1, psi1 = diagonalize_pmbdmatrix_gamma(np.conj(pqtoada).T @ mbdmatrix @ pqtoada, omegas)
        #print('eigs',eigs) 
        pt.exc += np.sum(eigs) / 4
        phi1 += phi1_t
        psi1 += psi1_t


    pt.exc = pt.exc / len(k_points) -3*np.sum(omega_rsscs) / 2 - np.sum(pt.omega) / 2
    phi1 = phi1 / len(k_points)
    psi1 = psi1 / len(k_points)
    pt.number, pt.adadaa, pt.X, pt.Y, pt.xx, pt.pp = calc_ptobservables(omegas, phi1, psi1, natm)

    if rpa is True:
      pt_AT_energy, pt_rpa_orders = get_AT_matrix(coords, alpha_0, C6, R_vdw, beta = beta, pt = pt, lattice=lattice, k_grid=k_grid, noint = noint, runmode = runmode, shift = shift)
    else:
      pt_AT_energy = 0
      pt_rpa_orders = 0

    print('ene per atoms: ', ene/natm)
    print('AT_energy per atom: ', AT_energy/natm)
    print('ptexc per atom: ', pt.exc/natm)
    print('pt_AT_energy per atom: ', pt_AT_energy/natm)

    output_array = {
       'ene' : ene,
       'ptexc' : pt.exc,
       'AT_energy' : AT_energy,
       'rpa_orders' : rpa_orders,
       'pt_AT_energy' : pt_AT_energy,
       'pt_rpa_orders' : pt_rpa_orders,
       'pt_number' : pt.number,
       'pt_adadaa' : pt.adadaa,
       'pt_X' : pt.X,
       'pt_Y' : pt.Y,
       'pt_xx' : pt.xx,
       'pt_pp' : pt.pp
    }

    return output_array

