import mass
import mass.materials
import numpy as np
import pylab as plt
from uncertainties import unumpy as unp
from uncertainties import ufloat 

xray_range = np.arange(500, 8000, 1)

EBIT_model = mass.materials.filterstack_models['EBIT 2018']
print(EBIT_model)
