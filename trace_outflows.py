import glob
import numpy as np 
import matplotlib
if __name__=="__main__":
	matplotlib.use('Agg') #makes sure matplotlib doesn't look for a display when plotting

import matplotlib.pylab as plt
import tangos
import pynbody 
from pynbody.plot import image
from matplotlib.pylab import cm
from scipy.stats import binned_statistic
import gc

cmap = cm.get_cmap('magma')
allgas = True
		
#load particles
basename = '/nobackupp2/mtremmel/Romulus/'
datadir = basename+'h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096'

#load halo database
db = tangos.get_simulation('h1.cosmo50')
steps = db.timesteps
stepname = [step.relative_filename for step in steps]
steptime = [step.time_gyr for step in steps]
stepnum = -1
halonum = 0 #generalize these later
halo = steps[stepnum].halos[halonum]
bh = halo['BH_central'][0]
halonum_pb = halo.halo_number
zs = np.array(halo.calculate_for_progenitors("z()")[0])
#total mass is in Msol

#average profiles
Tmw=halo.calculate_for_progenitors('Tmw_tcut_profile')[0]*pynbody.units.k.in_units("keV K**-1") ##converting k_B from J/K to keV/K, m_p from kg to g
rho_e_mw=halo.calculate_for_progenitors('rho_e_vol_profile')[0]
entropy_mw = Tmw*(rho_e_mw**(-2./3))
profiletime = halo.calculate_for_progenitors('t()')[0]

def entropy(ptcls):
	T_kev = ptcls['temp'].in_units('K')*pynbody.units.k.in_units('keV K**-1')
	if allgas:
		n_cc = ptcls['ne']*ptcls.g['rho'].in_units('m_p cm**-3')
	else:
		n_cc = ptcls.g['rho'].in_units('m_p cm**-3')/ptcls.g['mu']
	entropy = T_kev*pow(n_cc,-2./3)
	entropy.units = 'keV cm^2'
	return entropy #*weight

def ind(z):
	return np.argmin(abs(zs - z))

#select particles within certain distance of halo center
def filter(snap, halo_center=[0,0,0], radius='1 kpc'):
	snap.physical_units()
	spfilter = pynbody.filt.Sphere(radius, halo_center) & (pynbody.filt.FamilyFilter(pynbody.family.gas))# | pynbody.filt.FamilyFilter(pynbody.family.star))
	sphere = snap[spfilter]
	return sphere.g['iord']

def plot(halo_ptcls, agn_ptcls, snapnum, xmin=0, xmax=np.log10(500),
	ymin=-4.5, ymax=4, vmin=-6.5, vmax=-2.5, nbins=50, nosubs=True):
	fig, ax = plt.subplots()
	#get time stamp
	snap_halo = steps[snapnum].halos[halonum]
	tsnap = snap_halo.timestep.time_gyr
	
	#make histogram
	heatmap, xedges, yedges = np.histogram2d(np.log10(abs(halo_ptcls.g['r'])), 
		np.log10(entropy(halo_ptcls.g)),range=[[xmin, xmax],[ymin, ymax]], bins=nbins, weights = halo_ptcls.g['mass'])

	X, Y = np.meshgrid(xedges, yedges)
	hist = np.log10(heatmap.T/sum(heatmap.T)) #dividing by sum normalizes probability at each r
	plt.pcolormesh(X, Y, hist, vmin=-6.5, vmax=-2.5, cmap='gray_r')
	print "colormesh done"
	
	snap_halo = heatmap = X = Y = None #clear cache
	del(halo_ptcls)
	gc.collect()
	
	plt.colorbar()
	
	#track AGN
	agn_r = np.log10(abs(agn_ptcls.g['r']))
	plt.scatter(agn_r,np.log10(entropy(agn_ptcls.g)), color='r')
	print "AGN scatter done"
	
	snap_entropy = entropy_mw[np.argmin(abs(np.array(profiletime)-tsnap))] #because profiles calculated backwards from z=0
	profile_mask=np.ma.masked_invalid(snap_entropy).mask
	profile_r = (np.arange(len(snap_entropy))+1)/10.
	profile_r = np.ma.masked_array(profile_r, profile_mask).compressed()
	snap_entropy = np.ma.masked_array(snap_entropy, profile_mask).compressed()
	smooth = binned_statistic(np.log10(profile_r), snap_entropy, range=(0,np.log10(500)), bins=100)
	plt.plot(smooth.bin_edges[0:-1], np.log10(smooth.statistic), c=cmap(0.25), lw=2)
	print "overplot done"
	
	return fig, ax

snaps = glob.glob(datadir.split('.004096')[0]+'*')
unique_snaps = []
for snap in snaps:
	try:
		snap.split('1536gst1bwK1BH.')[1].split('.')[1]
	except IndexError:
		try: 
			float(snap.split('1536gst1bwK1BH.')[1].split('.')[0])
			unique_snaps.append(snap)
		except ValueError:
			continue

unique_snaps.sort()
peak_bh_activity = np.argmin(abs(np.array(steptime) - 6))

#bh_activity()


active_snap = pynbody.load(unique_snaps[peak_bh_activity])
halo_ptcls = active_snap.halos(dosort = True).load_copy(1)
halo_ptcls.physical_units()
pynbody.analysis.halo.center(halo_ptcls.g, mode='ssc')
cumul_indices = filter(halo_ptcls)
print "indices collected"
del(active_snap)# = None #clear cache
gc.collect()


def zoom(halo_ptcls, snapnum, halonum, snap, nosubs, cumul_indices=cumul_indices):
	new_indices = filter(halo_ptcls)
	cumul_indices = np.unique(np.concatenate([cumul_indices, new_indices]))
	agn_cumul_ptcls = halo_ptcls.g[(np.in1d(halo_ptcls.g['iord'], cumul_indices))]
	agn_cumul_ptcls.physical_units()

	fig, ax = plot(halo_ptcls, agn_cumul_ptcls, snapnum, xmin = 0, xmax = np.log10(2000), ymin = -2, 
		ymax = 3, nbins=50, nosubs=nosubs)

	suffix=''

	if nosubs:
		suffix += '_nosubs'

	xticks = np.array([5e1, 1e2, 5e2, 1e3])
	plt.xticks(np.log10(xticks), xticks)
	plt.xlabel("R (kpc)")
	plt.ylabel(r"K (KeV cm$^2$)")

	yticks = np.array([0.01, 0.1, 1, 10, 1e2, 1e3])
	plt.yticks(np.log10(yticks), yticks)
	plt.xlim(xmin, np.log10(2000))
	plt.ylim(-2, np.log10(500))
	plt.text(0.2, np.log10(2000),'%0.2f Gyr' % tsnap)
	plt.savefig('halo_%d_snap_%d_rain%s.png' % (halonum, snapnum, suffix))

	print "zoom done"

	del(agn_cumul_ptcls)
	#del(subhalo_ptcls)# = None
	gc.collect()


#if __name__=="__main__":
def trace_outflows(nosubs=True):
	#for snap in unique_snaps[57:]:
	for snapnum in [28, 41, 54,57, 60, 68]:#keysnaps: 
		#snapnum = unique_snaps.index(snap)
		snap = unique_snaps[snapnum]
		snap_halo = steps[snapnum].halos[halonum]
		stime = snap_halo.timestep.time_gyr
		print "t = %0.2f Gyr" % stime
		active_snap = pynbody.load(snap)
		halo_ptcls = active_snap.halos(dosort=True).load_copy(halonum_pb)
		if nosubs:
			halo_ptcls = halo_ptcls[halo_ptcls['amiga.grp'] == 1]
		halo_ptcls.physical_units()
		print "halo loaded"
		halo_ptcls['pos'] -= snap_halo['shrink_center']
		# # #2: track growing wind 
		zoom(halo_ptcls, snapnum, halonum, snap, nosubs, cumul_indices=cumul_indices)

		print "snap %d done" % snapnum
		del(active_snap)
		del(halo_ptcls)
		gc.collect()

