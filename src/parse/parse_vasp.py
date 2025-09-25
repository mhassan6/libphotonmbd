import re
import json
import numpy as np
from atomic_constants import P_Har, P_a_B

def parse_vasp_input(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    atomic_positions = []
    atomic_species = None
    lattice = []
    natom = None

    for cell_parameter_line in lines[2:5]:
        match = re.search(r'(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)', cell_parameter_line)
        if match:
            lattice.append([float(match.group(1)), float(match.group(2)), float(match.group(3))])

    atomic_species = re.findall(r'[A-Za-z]+', lines[5])

    natom = re.findall(r'(\d+)', lines[6])

    for i, line in enumerate(lines[8:]):
        match = re.search(r'(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)', line)
        if match:
            atomic_positions.append([float(match.group(1)), float(match.group(2)), float(match.group(3))])

    atomic_symbols = []
    for species, number in zip(atomic_species, natom):
        atomic_symbols += [species]*int(number)

    return {'atomic_positions':atomic_positions, 'atomic_species':atomic_symbols, 'lattice':lattice, 'natom':natom}

def parse_vasp_output(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    energy = []
    ratio = []
    natoms = None
    mbd_energy = None
    k_points = None
    xc = None

    for i, line in enumerate(lines):
        if "free energy" in line and "TOTEN" in line:
            match = re.search(r'free energy\s+TOTEN\s+=\s+(-?\d+\.\d+)\s+eV', line)
            if match:
                energy.append(float(match.group(1)))

    Edisp_flag = False
    for j, line in enumerate(lines):
        if "RELVOL" in line:
            for ratio_line in lines[j:]:
                match = re.search(r'[A-Za-z]+\s+(-?\d+\.\d+)', ratio_line)
                if match:
                    ratio.append(float(match.group(1)))

                if 'Edisp' in ratio_line:
                    Edisp_flag = True
                    break
        
        if Edisp_flag:
            break

    for line in lines:
        if "Edisp" in line:
            match = re.search(r'Edisp \(eV\):\s*(-?\d+\.\d+)', line)
            if match:
                mbd_energy = float(match.group(1))

    for line in lines:
        if "LIBMBD_XC" in line:
            match = re.search(r'LIBMBD_XC\s*=\s*(\w+)', line)
            if match:
                xc = match.group(1)
                
    for line in lines:
        if "generate k-points" in line:
            match = re.search(r'generate k-points for:\s+(\d+)\s+(\d+)\s+(\d+)', line)
            if match:
                k_points = [int(match.group(1)), int(match.group(2)), int(match.group(3))]

    return energy[-1], ratio, mbd_energy, k_points, xc

def parse_vasp(input_file, output_file, save_filename = 'vasp_input_output'):
    input_data = parse_vasp_input(input_file)
    energy, ratio, mbd_energy, k_points, xc = parse_vasp_output(output_file)
    
    results = {
    "code": 'VASP',
    "xc": xc,
    "atomic_species": input_data['atomic_species'],
    "atomic_positions": input_data['atomic_positions'],
    "lattice": input_data['lattice'],
    "kpoints": k_points,
    "total_energy": energy/P_Har,
    "ratio": ratio,
    "mbd_energy": mbd_energy/P_Har
    }

    with open(f"{save_filename}.json", "w") as f:
        json.dump(results, f, indent=2)