# libphotonmbd

libphotonmbd is a flexible software package that calculates cavity-modified many-body dispersion effects based on the Hirshfeld decomposition from a DFT calculation. As such, this package accepts DFT output files from different DFT codes, such as Quantum Espresso, VASP, and PySCF. Users will supply their DFT input and output files, which are then converted into a JSON file that is used internally by libphotonmbd to calculate photon-mbd energies. Regardless of the DFT code used, the resulting JSON file will always be formatted according to the standard found in the JSON.md markdown file.

