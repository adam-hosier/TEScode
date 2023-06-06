import mass
import mass.materials
import numpy as np
import pylab as plt
from uncertainties import unumpy as unp
from uncertainties import ufloat
import pandas as pd
dest = 'C:\\data'
xray_range = np.arange(200, 10000, 1)
EBIT_model_old = mass.materials.filterstack_models['EBIT 2018']
EBIT_model_new = mass.materials.FilterStack(name='2022 Stack')
EBIT_model_new.add_Film(name='Electroplated Au Absorber', material='Au',area_density_g_per_cm2=ufloat(0.00186,0.00006), fill_fraction=ufloat(1.000,0), absorber=True)
EBIT_model_new.add_AlFilmWithOxide(name='50mK Filter', Al_thickness_nm=100)
EBIT_model_new.add_AlFilmWithOxide(name='5K Filter', Al_thickness_nm=112)
EBIT_model_new.add_AlFilmWithOxide(name='50K Filter', Al_thickness_nm=116)
#EBIT_model_new.add_Film(name='Ni Mesh', material='Ni', area_density_g_per_cm2=ufloat(0.0134,0.0018), fill_fraction=ufloat(0.170,0.010), absorber=False)
EBIT_model_new.add_LEX_HT(name='Luxel Window #1')

eff = EBIT_model_new.get_efficiency(xray_range) 
eff_unc = EBIT_model_new.get_efficiency(xray_range, uncertain=True)
uncs = unp.std_devs(eff_unc)
# EBIT_model_old.plot_efficiency(xray_energies_eV=xray_range)
# plt.show()
#
#effdat = np.vstack((xray_range, eff, uncs)).T
effdat = np.vstack((xray_range, eff))
labels = ['Energy [eV]', 'Efficiency %']
for sur in EBIT_model_new.components.keys():
    effdat = np.vstack((effdat, EBIT_model_new.components.get(str(sur)).get_efficiency(xray_range)))
    labels.append(str(sur))

effdat = effdat.T 

EBIT_model_new.plot_efficiency(xray_energies_eV=xray_range)
plt.show()
#EBIT_model_old.plot_efficiency(xray_energies_eV=xray_range)

df = pd.DataFrame(data=effdat, columns=labels)
df.to_csv(dest+'/'+'TES_Efficiency_Dec2022.csv', index=False)
