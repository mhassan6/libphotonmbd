# Standard for JSON file

{  
  "code": "The DFT Code that was used",  
  "xc": "The exchange-correlation functional that was used",  
  "atomic_species": ["A list of the atomic species. For each atom, the atomic symbol is listed as many times as that atom appears in the structure"],  
  "atomic_positions": ["A list of the atomic position of each atom in Angstroms. Each position is a list that contains the x, y, and z positions"],  
  "lattice": ["The lattice dimensions in Angstroms"], 0 if the calculation was not periodic,  
  "kpoints": ["The k-grid given as a list of number of k-points in each direction"], 0 if the calculation was not periodic,  
  "total_energy": "The DFT energy, in Hartrees, from the DFT output file (excluding the mbd energy)",  
  "ratio": ["A list of the Hirshfeld ratio for each atom"],  
  "mbd_energy": "The mbd energy calculated through the DFT code, in Hartrees"  
}  
  
To illustrate this format, an example JSON file is provided for a VASP DFT calculation on a bilayer graphene unitcell:  
{  
  "code": "VASP",  
  "xc": "pbe",  
  "atomic_species": [  
    "C",  
    "C",  
    "C",  
    "C"  
  ],  
  "atomic_positions": [  
    [  
      0.0,  
      -0.0,  
      16.75000000000007  
    ],  
    [  
      1.2291642992972223,  
      0.710623905734562,  
      16.75000000000007  
    ],  
    [  
      0.0,  
      -0.0,  
      13.25000000000007  
    ],  
    [  
      1.2291642992972223,  
      0.710623905734562,  
      13.250000000000073  
    ]  
  ],  
  "lattice": [  
    [  
      2.46,  
      0.0,  
      0.0  
    ],  
    [  
      -1.23,  
      2.130422493309719,  
      0.0  
    ],  
    [  
      0.0,  
      -0.0,  
      30.0  
    ]  
  ],  
  "kpoints": [  
    1,  
    1,  
    1  
  ],  
  "total_energy": -0.8498713835717431,  
  "ratio": [  
    0.86744422,  
    0.86744745,  
    0.86744422,  
    0.86744745  
  ],  
  "mbd_energy": -0.10992678708990877  
}
