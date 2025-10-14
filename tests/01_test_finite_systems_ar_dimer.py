
import sys
sys.path.append("../src")
import numpy as np

from pypmbd_parse_json import pypmbd_parse_json
from photon_mode import photon_mode

from pytest import approx

pt = photon_mode(1, 0.2, 0.1, [0,0,1])
filename = 'json/ar_pyscf.json'
gauge = 'length'
output_array = pypmbd_parse_json(filename, pt, gauge)

vdwene_len = output_array['ene']
E0 = output_array['E0']
AT_energy = output_array['AT_energy']
rpa_orders = output_array['rpa_orders']
pt_AT_energy = output_array['pt_AT_energy']
pt_rpa_orders = output_array['pt_rpa_orders']
pt_number = output_array['pt_number']
pt_adadaa = output_array['pt_adadaa']

print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)

print('ene', vdwene_len, 'AT ene', AT_energy, 'ratio', vdwene_len/AT_energy)
print('log AT ene', AT_energy, 'rpa_order', np.sum(rpa_orders))

print('\n')
print('with photons')
print('ene', pt.exc, 'AT ene', pt_AT_energy, 'ratio', pt.exc/pt_AT_energy)
print('log AT ene', pt_AT_energy, 'rpa_order', np.sum(pt_rpa_orders))
print('photonnumber', pt_number)
print('photon adadaa', pt_adadaa)

#Compare to old code
#Gauge length E0 -1054.6393919338504 vdwene -0.00024635007299611544 vdw+pt 0.02928918951142423 pt 0.029535539584420345 ptnumber 0.010602342808323782 adadaa 0.001991878451505221 Mandel Q 0.177269195352225
print(approx(-0.00024635007299611544, abs=1e-10))
assert vdwene_len == approx(-0.00024635007299611544, abs=1e-10)
assert pt.exc == approx(0.02928918951142423, abs=1e-10)
assert pt_number == approx(0.010602342808323782, abs=1e-10)
assert pt_adadaa == approx(0.001991878451505221, abs=1e-10)

assert AT_energy == approx(-0.00024635007299611544, abs=1e-10)
assert pt_AT_energy == approx(0.02928918951142423, abs=1e-10)

print()
print('Length gauge tests passed!')
print()


gauge = 'velocity'
output_array = pypmbd_parse_json(filename, pt, gauge)

vdwene_len = output_array['ene']
E0 = output_array['E0']
AT_energy = output_array['AT_energy']
rpa_orders = output_array['rpa_orders']
pt_AT_energy = output_array['pt_AT_energy']
pt_rpa_orders = output_array['pt_rpa_orders']
pt_number = output_array['pt_number']
pt_adadaa = output_array['pt_adadaa']

print('Gauge', gauge, 'E0', E0, 'vdwene', vdwene_len, 'vdw+pt', pt.exc, 'pt', pt.exc-vdwene_len)

print('ene', vdwene_len, 'AT ene', AT_energy, 'ratio', vdwene_len/AT_energy)
print('log AT ene', AT_energy, 'rpa_order', np.sum(rpa_orders))

print('\n')
print('with photons')
print('ene', pt.exc, 'AT ene', pt_AT_energy, 'ratio', pt.exc/pt_AT_energy)
print('log AT ene', pt_AT_energy, 'rpa_order', np.sum(pt_rpa_orders))
print('photonnumber', pt_number)
print('photon adadaa', pt_adadaa)

#Gauge velocity E0 -1054.6393919338504 vdwene -0.00024635007299611544 vdw+pt 0.029289189511424674 pt 0.02953553958442079 ptnumber 0.11160152906407572 adadaa 0.04536898212856857 Mandel Q 0.2949249989239063
assert vdwene_len == approx(-0.00024635007299611544, abs=1e-10)
assert pt.exc == approx(0.02928918951142423, abs=1e-10)
assert pt_number == approx(0.11160152906407572, abs=1e-8)
assert pt_adadaa == approx(0.04536898212856857, abs=1e-10)

assert AT_energy == approx(-0.00024635007299611544, abs=1e-10)
assert pt_AT_energy == approx(0.02928918951142423, abs=1e-10)


print()
print('Velocity gauge tests passed!')
print()

#Reference data from pyqedft code used in C. Tasci et al., Phys. Rev. Lett. 134, 073002 (2025)
#mol.basis = 'ccpvdz'

#functional: pbe mbd_beta: 0.83
#converged SCF energy = -1054.63939193385
#This is not a periodic calculation
#Veff/Vfree ratio [1.00028763 1.00028763]
#ene -0.00024635007299611544 AT ene -0.0002463500726574358 ratio 1.0000000013747903
#log AT ene -0.0002463500726574358 rpa_order -0.000246350072657402
#ene 0.02928918951142423 AT ene 0.029289344157399172 ratio 0.9999947200601655
#log AT ene 0.029289344157399172 rpa_order 0.029289344157811686
#Gauge length E0 -1054.6393919338504 vdwene -0.00024635007299611544 vdw+pt 0.02928918951142423 pt 0.029535539584420345 ptnumber 0.010602342808323782 adadaa 0.001991878451505221 Mandel Q 0.177269195352225
#ene -0.00024635007299611544 AT ene -0.0002463500726574358 ratio 1.0000000013747903
#log AT ene -0.0002463500726574358 rpa_order -0.000246350072657402
#ene 0.029289189511424674 AT ene 0.029289344157399172 ratio 0.9999947200601806
#log AT ene 0.029289344157399172 rpa_order 0.029289344157811686
#Gauge velocity E0 -1054.6393919338504 vdwene -0.00024635007299611544 vdw+pt 0.029289189511424674 pt 0.02953553958442079 ptnumber 0.11160152906407572 adadaa 0.04536898212856857 Mandel Q 0.2949249989239063
