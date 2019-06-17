
#this analysis is officially easier to run on Grace
#because both the snapshots and database files now live here
#and pynbody and tangos work
import tangos as db
import numpy as np
import matplotlib.pylab as plt
from matplotlib import cm
import os
from scipy.stats import binned_statistic
from astropy import units, constants
from tangos.examples import mergers
# from astropy.stats import LombScargle

omega_m = 0.3
simname='h1.cosmo50'
sim = db.get_simulation(simname)
steps = sim.timesteps
#black holes count as halos
#get local copies of h1.cosmo50 and the database to the office desktop
haloind = 0
startstep = -1
h = steps[startstep].halos[0]
zs = np.array(h.calculate_for_progenitors("z()")[0])
ts = np.array(h.calculate_for_progenitors("t()")[0])
#total mass is in Msol
ind00 = np.argmin(abs(ts - 13.80)) #sim end
ind01 = np.argmin(abs(ts - 12.50)) #disruption end
ind02 = np.argmin(abs(ts - 12)) # disruption start
ind03 = np.argmin(abs(ts - 11.1)) #1st pericenter passage done
ind04 = np.argmin(abs(ts - 10.8)) #merger start
ind05 = np.argmin(abs(ts - 9.7))#post quenching end
ind06 = np.argmin(abs(ts - 9.2)) #post quenching start
ind07 = np.argmin(abs(ts - 7.7)) #pre-quenching end
ind08 = np.argmin(abs(ts - 7.3)) #pre-quenching start

epochs = list(range(ind01, ind02+1)) + list(range(ind03, ind04+1)) + list(range(ind05, ind06+1)) + list(range(ind07, ind08+1))
pairs = [[ind07,ind08], [ind05,ind06], [ind03, ind04], [ind01, ind02]]

cmap = cm.RdYlBu_r
def color(ind):
	if ind in range(ind07,ind08):           
		return cmap(0)
	elif ind in range(ind05, ind06):
		return cmap(0.2)
	elif ind in range(ind03, ind04):
		return cmap(0.75)
	elif ind in range(ind01, ind02):
		return cmap(0.95) 
	else:
		return cmap(0)

cummass_profile = h.reverse_property_cascade("tot_mass_profile")[0]
labels = {ind01:'Non-cool core', ind03:'Merger begins', ind05:'Isolated, loud AGN', ind07: 'Isolated, quiet AGN'}

rho = h.reverse_property_cascade('gas_density_profile')[0]

fM = 1.07 #McD 14 table 5
rho_crit=124.5786101694118#Msol/kpc^3
m_p = 1.67e-24#g
f_b = 0.14
mu = 0.59
A=1.397
Z=1.199
mu_e = A/Z

h70 = 0.67/.7
rho_crit_a = rho_crit*2e33/((3.085e21)**3)*(1 + zs)**3 #g/cm^3
ng_500 = 500*rho_crit_a*f_b/(mu*m_p) #mconald 2014 pg 7
ne_500 = ng_500*mu/mu_e

def r500_ind(row, overdensity=500, type='crit'):
	rho500 = overdensity*rho_crit_a[row]/(2e33/(3.086e21)**3)
	if type == 'mean':
		rho500 *= omega_m
	rsnap  = np.arange(len(cummass_profile[row]))/10.
	cumrho = np.divide(cummass_profile[row], 4*np.pi*(rsnap**3)/3.) #msol kpc^-3
	return np.argmin(abs(cumrho-rho500))

kT500s = (h.reverse_property_cascade('T500mw')[0] + h.reverse_property_cascade('T500ew')[0])/2.
P500 = ng_500*kT500s#mcdonald 2014 eq 4
K500 = kT500s*pow(ne_500,-2./3) #mcdonald 2014 eq 4
Tmw = h.reverse_property_cascade('Tmw_tcut_profile')[0]
Tmw *= constants.k_B.to('keV K**-1').value #for entropy to be in keV cm^2

rho_e_vol = np.array([rho*units.g/(units.cm**3) for rho in h.reverse_property_cascade('rho_e_vol_profile')[0]])
rho_g_vol = m_p*rho_e_vol*mu_e #g/cm^3

def entropy(T, rho_e):
	return T/pow(rho_e, 2/3.) 

Kmw = np.array([entropy(Tmw[elt].in_units('keV')*units.keV, rho_g_vol[elt]/(m_p*units.g)) for elt in range(len(Tmw))])
profilelist = {r'K/K$_{500}$':Kmw, r'$\rho/\rho_{crit}$':rho_g_vol, r'T/T$_{500}$':Tmw}#, 'precip':-precip}
norms = {r'K/K$_{500}$':K500, r'$\rho/\rho_{crit}$':rho_crit_a, r'T/T$_{500}$':kT500s}#, 'precip':None}
T500s = h.reverse_property_cascade('T500mw')[0]*units.keV
cs = np.sqrt((T500s/(mu*constants.m_p)).to('kpc**2 yr**-2'))

bh = h['BH_central'][0]
bh_mdot = bh.reverse_property_cascade('BH_mdot_ave')[0]
bh_ts = bh.calculate_for_progenitors('t()')[0]
bh_mass = bh.reverse_property_cascade('BH_mass')[0]

def Pbh():
	eps = 0.02
	#sound crossing time - no I really need adiabatic index for this
	power = eps*(bh_mdot)*units.Msun*constants.c**2/units.yr
	return power.value*power.unit.in_units('keV yr**-1')

def r_E(ind): #a la Tang & Churazov
	t = tbh[ind]
	mdot = bhmdot[ind]
	temp = Tmw/pynbody.units.k.in_units('keV K**-1') #back to K
	rho_e = h['rho_e_vol_profile'] #cm^-3
	rho_g = rho_e/mu
	Etherm_vol = 3./2 * pynbody.units.k.in_units('erg K**-1')*rho_g*units.cc.in_units('m^-3')*temp

	dr = 0.1*pynbody.units.kpc.in_units('m')#kpc
	r = np.arange(len(rho_e))/10. *pynbody.units.kpc.in_units('m') #kpc 
	dEtherm = Etherm_vol * 4 *np.pi*r**2 * dr 
	Etherm = np.cumsum(dEtherm) #units of ergs

def Pnorm(x, P0, c500, alpha, gamma): #Arnaud et al 2010 M500
    beta = 5.4905
    ap = 0.12 #eq 7
    apx = 0.1 - (ap+0.10)*pow(x/0.5,3)/(1+pow(x/0.5,3)) #eq 8
    px = P0/(pow(c500*x, gamma)*pow(1+pow(c500*x, alpha), (beta-gamma)/alpha)) #eq 11
    return px#*pow(M500*h70/3e14,0.12) #eq 13

def bin_data(profile, elt, overdensity = 500, type='crit', ynorm=1, xnorm=False,statistic='mean'):
	if xnorm:
		r500 = r500_ind(elt, overdensity, type)
		print( r500/10., 'kpc = R500')
		bins=np.linspace(np.log10(.003), np.log10(2), 101)
	else:
		r500 = 10.
		bins=np.linspace(0, 3, 101)
	profile_mask = np.ma.masked_invalid(profile[elt]).mask
	rsnap = np.arange(len(profile[elt]))
	profile_r = np.ma.masked_array(rsnap, profile_mask).compressed()/float(r500)
	profile_r+=0.001
	snap_profile = np.ma.masked_array(profile[elt], profile_mask).compressed()/ynorm[elt]
	smooth = binned_statistic(np.log10(profile_r), snap_profile, bins=bins, statistic=statistic)
	return smooth

def profile_array(profile, ind1, ind2, overdensity=500, type='crit',ynorm=None, xnorm=True, statistic='mean'):
	indices = range(ind1, ind2)
	profile_array = np.empty([len(indices), 100])
	for elt in indices:
		smooth = bin_data(profile, elt, overdensity, type, ynorm=ynorm, xnorm=xnorm, statistic=statistic)
		profile_array[elt-ind1] = smooth.statistic
	return profile_array 

def plot_profile(profile_array, color, norm=True, label=None, fill=True, linestyle='solid', ax=None):
    if norm:
    	bins = np.linspace(np.log10(.003), np.log10(2), 100)
    else:
    	bins = np.linspace(0, 3, 100)
    mean = np.nanmean(profile_array,axis=0)
    std = np.nanstd(profile_array,axis=0)
    min = np.nanmin(profile_array,axis=0)
    if ax == None:
    	fig, ax = plt.subplots()
    	ax.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
    	ax.xaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
    if fill:
	    ax.fill_between(10**bins, np.maximum(mean-std,min), mean+std, color=color,alpha=0.5)
    ax.plot(10**bins, mean, color=color, label=label,linestyle=linestyle)