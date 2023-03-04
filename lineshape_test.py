import mass 
import matplotlib.pylab as plt 
import numpy as np 

line = mass.FeKBeta
#line2 = mass.SKBeta 

N = 100000
energies = line.rvs(size=N, instrument_gaussian_fwhm=0)
plt.clf() 
sim, bin_edges, _ = plt.hist(energies, 300, [7000, 7150], histtype="step");
binsize = bin_edges[1] - bin_edges[0]
e = bin_edges[:-1] + 0.5*binsize
plt.plot(e, line(e, instrument_gaussian_fwhm=2.2)*N*binsize, "k")
plt.xlabel("Energy (eV)")
plt.show() 