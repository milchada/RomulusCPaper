
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
