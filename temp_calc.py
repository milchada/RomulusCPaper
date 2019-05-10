from astropy import constants, units
import pynbody 
import numpy as np

def emissivity(rho, T, mu, tcool):
	#taken from https://arxiv.org/abs/astro-ph/9809159 eq 11
	#mu = 0.59 #assumed mean molecular weight
	return 3./2. * rho * constants.k_B.value * T  / (mu*constants.m_p.value*tcool) #1e4 to convert m^2 in kB to cm^2
	#units shouldn't matter because we normalize over sum

def tcool(rho, T, mu):
	#taken from https://arxiv.org/abs/astro-ph/9809159 eq 12
	fm = 1.0 #metalicity dependent factor, 1.0 for solar, 0.03 for pristine
	# mu = 0.59
	C1 = 3.88e11
	C2 = 5e7
	return C1*mu*constants.m_p.value * T**0.5/(rho*(1+C2*fm/T)) #1e3 to convert m_p from kg to g 
	#units of s

def cluster_temp(h,r,nocore=False, ew=True, temp_cut=None):
	if temp_cut is None:
		temp_cut = 0
	if nocore is True:
		#this case is buggy
		use = np.where((h.g['temp'] > temp_cut) & (h.g['r'] < r)  & (h.g['r'] > 0.15 * r))[0]
	else:
		use = np.where((h.g['temp'] > temp_cut) & (h.g['r'] < r) )[0]

	if ew:
		mu = 0.58
		tc = tcool(h.g['rho'][use].in_units('g cm**-3'),h.g['temp'][use], mu)
		em = emissivity(h.g['rho'][use].in_units('g cm**-3'),h.g['temp'][use],mu,tc)
		Tew = np.sum(em*h.g['temp'][use])/np.sum(em)
		return constants.k_B.value*units.J.to('keV')*Tew #converting k_B from J/K to keV/K

	else:
		Tmw = np.sum(h.g['mass'][use]*h.g['temp'][use])/np.sum(h.g['mass'][use])
		return constants.k_B.value*units.J.to('keV')*Tmw

	#returns temp in keV

# #for some reason i can't read gas properties for the central halos
def main_cluster_temp(nocore=True, ew=True, den_cut=None,temp_cut=None):
	h1=s.halos(dosort=True)[1]
	r =h1.properties.get('Rvir')

	if den_cut is None:
		rho_cut = h1.properties.get('rho').max()
	else:
		rho_c = pynbody.analysis.cosmology.rho_crit(h, z=0) * (1.0 + h1.properties.get("z")) ** 3
		rho_cut = den_cut * rho_c
	if temp_cut is None:
		temp_cut = 0
	if nocore is True:
		use = np.where((h1.properties.get('temp') > temp_cut) & (h1.properties.get('r') < r) & (h1.properties.get('rho') <= rho_cut) & (h1.properties.get('r') > 0.15 * r))[0]
	else:
		use = np.where((h1.properties.get('temp') > temp_cut) & (h1.properties.get('r') < r) & (h1.properties.get('rho') <= rho_cut))[0]

	if ew:
		mu = 0.58
		tc = tcool(h1.properties.get('rho')[use].in_units('g cm**-3'),h1.properties.get('temp')[use], mu)
		em = emissivity(h1.properties.get('rho')[use].in_units('g cm**-3'),h1.properties.get('temp')[use],mu,tc)
		Tew = np.sum(em*h['temp'][use])/np.sum(em)
		return constants.k_B*units.J.to('keV')*Tew #converting k_B from J/K to keV/K

	else:
		Tmw = np.sum(h1.properties.get('mass')[use]*h1.properties.get('temp')[use])/np.sum(h1.properties.get('mass')[use])
		return constants.k_B*units.J.to('keV')*Tmw