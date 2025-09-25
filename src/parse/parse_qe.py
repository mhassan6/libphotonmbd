import re
import json
import numpy as np

def parse_qe_output(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    vfree = []
    veff = []
    natoms = None
    mbd_energy = []

    for i, line in enumerate(lines):
        if "!" in line and "total energy" in line:
            match = re.search(r'total energy\s+=\s+(-?\d+\.\d+)', line)
            if match:
                energy = float(match.group(1))

    for line in lines:
        if "vfree_pub" in line:
            match = re.search(r'vfree_pub for atom #\s*\d+\s*is\s*(-?\d+\.\d+)\s*bohr\^3 and veff_pub is\s*(-?\d+\.\d+)\s*bohr\^3', line)
            if match:
                vfree.append(float(match.group(1)))
                veff.append(float(match.group(2)))

    for line in lines:
        if "number of atoms/cell" in line:
            match = re.search(r'number of atoms/cell\s*=\s*(\d+)', line)
            if match:
                natoms = float(match.group(1))


    for line in lines:
        if "MBD Energy (EmbdvdW)" in line:
            match = re.search(r'MBD Energy \(EmbdvdW\) is\s*(-?\d+\.\d+)\s*Ry.', line)
            if match:
                mbd_energy.append(float(match.group(1)))

    veff = veff[len(veff)-int(natoms):]

    return energy, vfree, veff, mbd_energy[-1]

def parse_qe_input(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    atomic_positions = []
    atomic_species = []
    lattice = []
    kpoints = None
    xc = None

    for i, line in enumerate(lines):
        if "ATOMIC_POSITIONS" in line:
            for atomic_position_line in lines[i+1:]:
                match = re.search(r'(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)', atomic_position_line)
                if match:
                    atomic_species.append(match.group(1))
                    atomic_positions.append([float(match.group(2)), float(match.group(3)), float(match.group(4))])

    for i, line in enumerate(lines):
        if "CELL_PARAMETERS" in line:
            for cell_parameter_line in lines[i+1:i+4]:
                match = re.search(r'(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)', cell_parameter_line)
                if match:
                    lattice.append([float(match.group(1)), float(match.group(2)), float(match.group(3))])

    for i, line in enumerate(lines):
        if "K_POINTS" in line:
            match_kpoints = re.search(r'(\d+)\s*(\d+)\s*(\d+)\s*(\d+)\s*(\d+)\s*(\d+)', lines[i+1])
            match_option = re.search(r'K_POINTS\s*(\w+)',line)
            if match_kpoints:
                kpoints = [int(match_kpoints.group(1)), int(match_kpoints.group(2)), int(match_kpoints.group(3))]
            elif match_option:
                if match_option.group(1) == 'gamma':
                    kpoints = [1,1,1]
    
    for line in lines:
        if "input_dft" in line:
            match = re.search(r'input_dft\s*=\s*\'(\w+)\'', line)
            if match:
                xc = match.group(1)


    return atomic_positions, atomic_species, lattice, kpoints, xc


def parse_qe(input_file, output_file, save_filename = 'qe_input_output'):
    atomic_positions, atomic_species, lattice, kpoints, xc = parse_qe_input(input_file)
    energy, vfree, veff, mbd_energy = parse_qe_output(output_file)
    ratio = np.array(veff)/np.array(vfree)

    results = {
    "code": 'QE',
    "xc": xc,
    "atomic_species": atomic_species,
    "atomic_positions": atomic_positions,
    "lattice": lattice,
    "kpoints": kpoints,
    "total_energy": 0.5*energy,
    "ratio": ratio,
    "mbd_energy": mbd_energy
    }

    with open(f"{save_filename}.json", "w") as f:
        json.dump(results, f, indent=2)

