#!/bin/bash

python3 01_test_finite_ar_dimer.py | tail -1
python3 02_test_finite_ar_dimer_two_modes.py | tail -1
python3 03_test_finite_ar_dimer_mbd_AT_rpa_orders.py  | tail -1
python3 04_periodic_ar_mbd.py  | tail -1
python3 05_periodic_ar_pt_supercell.py | tail -1
python3 06_periodic_ar_pt_unitcell.py  | tail -1
