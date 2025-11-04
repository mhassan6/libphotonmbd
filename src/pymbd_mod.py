# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as np

from pymbd.pymbd import freq_grid, dipole_matrix, _array, screening, get_kpts
from pymbd.fortran import MBDGeom

#
# Adapted screening function from pymbd
# (https://github.com/libmbd/libmbd/tree/master/src/pymbd)
#

def screening_all_freqs(coords, alpha_0, C6, R_vdw, beta, lattice=None, nfreq=30):
    r"""Screen atomic polarizabilities.

    :param array-like coords: (a.u.) atom coordinates in rows
    :param array-like alpha_0: (a.u.) atomic polarizabilities
    :param array-like C6: (a.u.) atomic :math:`C_6` coefficients
    :param array-like R_vdw: (a.u.) atomic vdW radii
    :param float beta: MBD damping parameter :math:`\beta`
    :param array-like lattice: (a.u.) lattice vectors in rows
    :param int nfreq: number of grid points for frequency quadrature

    Returns static polarizabilities, :math:`C_6` coefficients, and
    :math:`R_\mathrm{vdw}` coefficients (a.u.).
    """
    freq, freq_w = freq_grid(nfreq)
    omega = 4 / 3 * C6 / alpha_0**2
    alpha_dyn = [alpha_0 / (1 + (u / omega) ** 2) for u in freq]
    alpha_dyn_rsscs = []
    alphanewr = []
    for a in alpha_dyn:
        sigma = (np.sqrt(2 / np.pi) * a / 3) ** (1 / 3)
        dipmat = dipole_matrix(
            coords, 'fermi,dip,gg', sigma=sigma, R_vdw=R_vdw, beta=beta, lattice=lattice
        ) #fermi,dip,gg
        a_nlc = np.linalg.inv(np.diag(np.repeat(1 / a, 3)) + dipmat)
        a_contr = sum(np.sum(a_nlc[i::3, i::3], 1) for i in range(3))/3
        a1=np.sum(a_nlc[0::3, 0::3], 1)
        a2=np.sum(a_nlc[1::3, 1::3], 1)
        a3=np.sum(a_nlc[2::3, 2::3], 1)
        a_contrnew = np.array([a1,a2,a3])
        alphanewr.append(a_contrnew)
        alpha_dyn_rsscs.append(a_contr)
    alphanewr=np.stack(alphanewr)
    alpha_dyn_rsscs = np.stack(alpha_dyn_rsscs)
    C6_rsscs = 3 / np.pi * np.sum(freq_w[:, None] * alpha_dyn_rsscs**2, 0)
    R_vdw_rsscs = R_vdw * np.sign(alpha_dyn_rsscs[0, :] / alpha_0)*np.abs(alpha_dyn_rsscs[0, :] / alpha_0)** (1 / 3)
    C6_rsscsnewr_a1 = 3 / np.pi * np.sum(freq_w[:, None] * alphanewr[:,0,:]**2, 0)
    C6_rsscsnewr_a2 = 3 / np.pi * np.sum(freq_w[:, None] * alphanewr[:,1,:]**2, 0)
    C6_rsscsnewr_a3 = 3 / np.pi * np.sum(freq_w[:, None] * alphanewr[:,2,:]**2, 0)
    C6_contrnew = np.array([C6_rsscsnewr_a1,C6_rsscsnewr_a2,C6_rsscsnewr_a3])
    omega_rsscs_a1 = 4 / 3 * C6_rsscsnewr_a1 / alphanewr[0,0,:]**2
    omega_rsscs_a2 = 4 / 3 * C6_rsscsnewr_a2 / alphanewr[0,1,:]**2
    omega_rsscs_a3 = 4 / 3 * C6_rsscsnewr_a3 / alphanewr[0,2,:]**2
    omega_contrnew = np.array([omega_rsscs_a1, omega_rsscs_a2, omega_rsscs_a3])
    return alpha_dyn_rsscs,alphanewr[0],alphanewr,freq,freq_w, omega_contrnew

# this is the same function as in pymbd -> md_energy(), only with output of all frequencies
def screening_all_freqs_orig(coords, alpha_0, C6, R_vdw, beta, lattice=None, nfreq=15):
    r"""Screen atomic polarizabilities.

    :param array-like coords: (a.u.) atom coordinates in rows
    :param array-like alpha_0: (a.u.) atomic polarizabilities
    :param array-like C6: (a.u.) atomic :math:`C_6` coefficients
    :param array-like R_vdw: (a.u.) atomic vdW radii
    :param float beta: MBD damping parameter :math:`\beta`
    :param array-like lattice: (a.u.) lattice vectors in rows
    :param int nfreq: number of grid points for frequency quadrature

    Returns static polarizabilities, :math:`C_6` coefficients, and
    :math:`R_\mathrm{vdw}` coefficients (a.u.).
    """
    freq, freq_w = freq_grid(nfreq)
    omega = 4 / 3 * C6 / alpha_0**2
    alpha_dyn = [alpha_0 / (1 + (u / omega) ** 2) for u in freq]
    alpha_dyn_rsscs = []
    for a in alpha_dyn:
        sigma = (np.sqrt(2 / np.pi) * a / 3) ** (1 / 3)
        dipmat = dipole_matrix(
            coords, 'fermi,dip,gg', sigma=sigma, R_vdw=R_vdw, beta=beta, lattice=lattice
        )
        a_nlc = np.linalg.inv(np.diag(np.repeat(1 / a, 3)) + dipmat)
        a_contr = sum(np.sum(a_nlc[i::3, i::3], 1) for i in range(3)) / 3
        alpha_dyn_rsscs.append(a_contr)
    alpha_dyn_rsscs = np.stack(alpha_dyn_rsscs)
    C6_rsscs = 3 / np.pi * np.sum(freq_w[:, None] * alpha_dyn_rsscs**2, 0)
    R_vdw_rsscs = R_vdw * (alpha_dyn_rsscs[0, :] / alpha_0) ** (1 / 3)
    return alpha_dyn_rsscs, C6_rsscs, R_vdw_rsscs

# this is the same function as in pymbd -> md_energy(), only with additional output
def mbd_energy_get_matrix(coords, alpha_0, C6, R_vdw, beta, lattice=None, k_grid=None, nfreq=15, noint = False, shift=0.5):
    r"""Calculate an MBD energy.

    :param array-like coords: (a.u.) atom coordinates in rows
    :param array-like alpha_0: (a.u.) atomic polarizabilities
    :param array-like C6: (a.u.) atomic :math:`C_6` coefficients
    :param array-like R_vdw: (a.u.) atomic vdW radii
    :param float beta: MBD damping parameter :math:`\beta`
    :param array-like lattice: (a.u.) lattice vectors in rows
    :param array-like k_grid: number of :math:`k`-points along reciprocal axes
    :param int nfreq: number of grid points for frequency quadrature
    """
    coords, alpha_0, C6, R_vdw, lattice = map(
        _array, (coords, alpha_0, C6, R_vdw, lattice)
    )
    alpha_0_rsscs, C6_rsscs, R_vdw_rsscs = screening(
        coords, alpha_0, C6, R_vdw, beta, lattice=lattice, nfreq=15
    )
    print('alpha_0_rsscs: ', alpha_0_rsscs)
    omega_rsscs = 4 / 3 * C6_rsscs / alpha_0_rsscs**2
    
    pre = np.repeat(omega_rsscs * np.sqrt(alpha_0_rsscs), 3)

    if lattice is None:
        k_points = [None]
    else:
        assert k_grid is not None
        k_points = get_kpts(lattice, k_grid, shift = shift)
    matrix_list = []
    ene = 0.0
    for k_point in k_points:
        matrix = np.diag(np.repeat(omega_rsscs**2, 3))

        dipole_matrix_mbd = dipole_matrix(
                coords,
                'fermi,dip',
                R_vdw=R_vdw_rsscs,
                beta=beta,
                lattice=lattice,
                k_point=k_point,
            )

        if noint:
          dipole_matrix_mbd *= 0
        
        matrix = np.diag(np.repeat(omega_rsscs**2, 3)) + np.outer(pre, pre) \
            * dipole_matrix_mbd
        
        matrix_list.append(matrix)

        # eigs = np.linalg.eigvalsh(matrix)
        eigs, eigvecs = np.linalg.eigh(matrix)

        
        ene += np.sum(np.sqrt(eigs)) / 2
    ene = ene / len(k_points) - 3 * np.sum(omega_rsscs) / 2

    return ene, alpha_0_rsscs, omega_rsscs, matrix_list#, dipole_matrix_mbd
        
# this is the same function as in pymbd -> md_energy(), only with additional output, but more performant

def mbd_energy_get_matrix2(coords, alpha_0, C6, R_vdw, beta, lattice=None, k_grid=None, nfreq=15, noint = False, runmode = 'rsscs', shift=0.5):
    r"""Calculate an MBD energy.

    :param array-like coords: (a.u.) atom coordinates in rows
    :param array-like alpha_0: (a.u.) atomic polarizabilities
    :param array-like C6: (a.u.) atomic :math:`C_6` coefficients
    :param array-like R_vdw: (a.u.) atomic vdW radii
    :param float beta: MBD damping parameter :math:`\beta`
    :param array-like lattice: (a.u.) lattice vectors in rows
    :param array-like k_grid: number of :math:`k`-points along reciprocal axes
    :param int nfreq: number of grid points for frequency quadrature
    """
    coords, alpha_0, C6, R_vdw, lattice = map(
        _array, (coords, alpha_0, C6, R_vdw, lattice)
    )
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

    pre = np.repeat(omega_rsscs * np.sqrt(alpha_0_rsscs), 3)
    pre2 = np.outer(pre, pre)
    if lattice is None:
        k_points = [None]
        Mgeom = MBDGeom(coords)
    else:
        assert k_grid is not None
        k_points = get_kpts(lattice, k_grid, shift = shift)
        Mgeom = MBDGeom(coords, lattice, k_points)        

    matrix_list = []
    eigvecs_list = []
    eigs_list = []

    ene = 0.0
    matrix_diag = np.diag(np.repeat(omega_rsscs**2, 3))

    if runmode == 'plain2':
      sigma_array = []
      for a in alpha_0_rsscs:
        sigma = (np.sqrt(2 / np.pi) * a / 3) ** (1 / 3)
        sigma_array.append(sigma)

    for k_point in k_points:
        if not noint:
          if runmode == 'plain2':
            dipole_matrix_mbd = Mgeom.dipole_matrix('dip,gg', sigma=sigma_array, k_point=k_point)
          else:
            dipole_matrix_mbd = Mgeom.dipole_matrix('fermi,dip', beta=beta, k_point=k_point, R_vdw= R_vdw_rsscs)

          #dipole_matrix_mbd = dipole_matrix(
          #        coords,
          #        'fermi,dip',
          #        R_vdw=R_vdw_rsscs,
          #        beta=beta,
          #        lattice=lattice,
          #        k_point=k_point)
        else:
          dipole_matrix_mbd = np.zeros((3*natm, 3*natm))

        matrix = matrix_diag + pre2 \
            * dipole_matrix_mbd

        #eigs = np.linalg.eigvalsh(matrix)
        eigs, eigvecs = np.linalg.eigh(matrix)
        ene += np.sum(np.sqrt(eigs)) / 2

        matrix_list.append(matrix)
        eigvecs_list.append(eigvecs)
        eigs_list.append(eigs)       
 
    ene = ene / len(k_points) - 3 * np.sum(omega_rsscs) / 2

    return ene, alpha_0_rsscs, omega_rsscs, matrix_list, eigvecs_list, eigs_list #, dipole_matrix_mbd
