import numpy as np
from read_db import Tmw, rho_e_vol, ts

def slope(pressure, density):
    pvalid = np.ma.masked_invalid(pressure).compressed()
    rhovalid = np.ma.masked_invalid(density).compressed()
    return np.polyfit(np.log10(rhovalid),np.log10(pvalid),1)

def polytropic_index(in_rkpc=len(Tmw[0])/10., label=None, filename=None):
	in_rbin = int(in_rkpc*10)
	pressure = Tmw*rho_e_vol
	gamma = np.array([slope(pressure[row][:in_rbin],rho_e_vol[row][:in_rbin])[0] for row in xrange(len(Tmw))])
	if label != None:
		plt.clf()
		plt.plot(ts, gamma, label=label)
		plt.xlabel('Time (Gyr)')
		plt.ylabel(r'$\gamma$')
		plt.savefig('polytropic_index'+filename+'.png')