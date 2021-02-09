import glob
import numpy as np 
import matplotlib
if __name__=="__main__":
	matplotlib.use('Agg') #makes sure matplotlib doesn't look for a display when plotting

import matplotlib.pylab as plt
import tangos
import pynbody 
from tangos.examples import mergers
from pynbody.plot import image
from matplotlib.pylab import cm
from scipy.stats import binned_statistic
import gc

plot_rain = False
plot_zoom = False
plot_mixing = True
cmap = cm.get_cmap('magma')
		
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
mtree = mergers.get_mergers_of_major_progenitor(halo)
h2=mtree[2][3][1]
merger_snap = basename+h2.path.split('/halo')[0]
merger_sim = pynbody.load(merger_snap)
# h2_ptcls = s.halos(dosort=True).load_copy(h2.halo_number)
# subhalo_indices = h2_ptcls.get_index_list(s)
zs = np.array(halo.calculate_for_progenitors("z()")[0])
#total mass is in Msol

if plot_zoom:
	Tmw=halo.calculate_for_progenitors('Tmw_tcut_profile')[0]*pynbody.units.k.in_units("keV K**-1") ##converting k_B from J/K to keV/K, m_p from kg to g
	rho_e_mw=halo.calculate_for_progenitors('rho_e_vol_profile')[0]
	entropy_mw = Tmw*(rho_e_mw**(-2./3))
	profiletime = halo.calculate_for_progenitors('t()')[0]

def entropy(ptcls, mw=False):
	T_kev = ptcls['temp'].in_units('K')*pynbody.units.k.in_units('keV K**-1')
	n_cc = ptcls['ne']*ptcls.g['rho'].in_units('m_p cm**-3')
	entropy = T_kev*pow(n_cc,-2./3)
	entropy.units = 'keV cm^2'
	return entropy #*weight

def ind(z):
	return np.argmin(abs(zs - z))

#select particles within certain distance of halo center
def filter(snap=merger_sim,halo_center=[0,0,0], radius='1 kpc'):
	snap.physical_units()
	spfilter = pynbody.filt.Sphere(radius, halo_center) & (pynbody.filt.FamilyFilter(pynbody.family.gas))# | pynbody.filt.FamilyFilter(pynbody.family.star))
	sphere = snap[spfilter]
	# indices = sphere.get_index_list(snap)
	return sphere.g['iord']

def plot(halo_ptcls, agn_ptcls, snapnum, suffix='', xmin=0, xmax=np.log10(500),
	ymin=-4.5, ymax=4, yspace=20, nbins=50, overplot=False, subhalo_ptcls=None, 
	track_agn = True, track_subhalo_ptcls = False):
	plt.clf()
	#get time stamp
	snap_halo = steps[snapnum].halos[halonum]
	tsnap = snap_halo.timestep.time_gyr
	#make histogram
	heatmap, xedges, yedges = np.histogram2d(np.log10(abs(halo_ptcls.g['r'])), 
		np.log10(entropy(halo_ptcls.g)),range=[[xmin, xmax],[ymin, ymax]], bins=nbins, weights = halo_ptcls.g['mass'])

	X, Y = np.meshgrid(xedges, yedges)
	hist = np.log10(heatmap.T/sum(sum(heatmap.T))) #dividing by sum normalizes probability at each r
	plt.pcolormesh(X, Y, hist, vmin=-5.5, vmax=-1.5, cmap='gray')
	print "colormesh done"
	snap_halo = heatmap = X = Y = None #clear cache
	del(halo_ptcls)
	gc.collect()
	
	plt.colorbar()
	if track_agn:
		agn_r = np.log10(abs(agn_ptcls.g['r']))
		plt.scatter(agn_r,np.log10(entropy(agn_ptcls.g)), color='r')
		print "AGN scatter done"
	if track_subhalo_ptcls:
		subhalo_r = np.log10(abs(subhalo_ptcls.g['r']))
		plt.scatter(subhalo_r, np.log10(entropy(subhalo_ptcls.g)), color='b', alpha=0.1)
		print "subhalo scatter done"
	xticks = np.array([1,10,1e2,5e2])
	plt.xticks(np.log10(xticks), xticks)
	plt.xlabel("R (kpc)")
	plt.ylabel(r"K (KeV cm$^2$)")
	if overplot:
		suffix+='_zoom'
		snap_entropy = entropy_mw[np.argmin(abs(np.array(profiletime)-tsnap))] #because profiles calculated backwards from z=0
		profile_mask=np.ma.masked_invalid(snap_entropy).mask
		profile_r = (np.arange(len(snap_entropy))+1)/10.
		profile_r = np.ma.masked_array(profile_r, profile_mask).compressed()
		snap_entropy = np.ma.masked_array(snap_entropy, profile_mask).compressed()
		smooth = binned_statistic(np.log10(profile_r), snap_entropy, range=(0,np.log10(500)), bins=100)
		plt.plot(smooth.bin_edges[0:-1], np.log10(smooth.statistic), c=cmap(0.25), lw=2)
		print "overplot done"

	yticks = np.array([0.01,0.1,1,10,1e2,1e3])
	plt.yticks(np.log10(yticks), yticks)
	plt.xlim(xmin, np.log10(5e2))
	plt.ylim(-2, np.log10(300))
	plt.text(0.2, np.log10(2000),'%0.2f Gyr' % tsnap)
	plt.savefig('halo_%d_snap_%d_rain%s.png' % (halonum, snapnum, suffix))

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
if plot_rain or plot_mixing:
	merger_sim = pynbody.load('/nobackupp2/mtremmel/Romulus/h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.003360')
	if plot_rain:
		subhalo_ptcls = merger_sim.halos(dosort=True).load_copy(2)
		subhalo_indices = subhalo_ptcls.g['iord']
	if plot_mixing:
		halo_ptcls = merger_sim.halos(dosort = True).load_copy(1)
		halo_ptcls.physical_units()
		pynbody.analysis.halo.center(halo_ptcls.g, mode='ssc')
		core_ind = filter(halo_ptcls, radius='30 kpc')

	del(merger_sim, halo_ptcls)# = None #clear cache

if plot_zoom:
    active_snap = pynbody.load(unique_snaps[peak_bh_activity])
    halo_ptcls = active_snap.halos(dosort = True).load_copy(1)
    halo_ptcls.physical_units()
    pynbody.analysis.halo.center(halo_ptcls.g, mode='ssc')
    cumul_indices = filter(halo_ptcls)
    print "indices collected"
    del(active_snap)# = None #clear cache
    gc.collect()

def rain(active_snap, halo_ptcls,  snapnum):
	subhalo_ptcls = active_snap.g[(np.in1d(active_snap.g['iord'], subhalo_indices))]
	subhalo_ptcls.physical_units()
	snap_halo = steps[snapnum].halos[halonum]
	#subhalo_ptcls['pos'] -= snap_halo['shrink_center']
	agn_ptcls = None
	plot(halo_ptcls, agn_ptcls,snapnum, ymin = -2, 
		ymax = max(np.log10(entropy(halo_ptcls.g))), nbins=50, subhalo_ptcls=subhalo_ptcls, 
		track_agn=False, track_subhalo_ptcls = True)
	print "rain plot done"
	del(subhalo_ptcls) # = None
	gc.collect()

def zoom(halo_ptcls, steps, snapnum, halonum, snap, cumul_indices):
	new_indices = filter(halo_ptcls)
	cumul_indices = np.unique(np.concatenate([cumul_indices, new_indices]))
	agn_cumul_ptcls = halo_ptcls.g[(np.in1d(halo_ptcls.g['iord'], cumul_indices))]
	agn_cumul_ptcls.physical_units()
	#subhalo_ptcls = active_snap.g[(np.in1d(active_snap.g['iord'], subhalo_indices))]
	#subhalo_ptcls.physical_units()
	# print "wind indices collected"
	plot(halo_ptcls, agn_cumul_ptcls, snapnum, xmin = 0, xmax = np.log10(2000), ymin = -2, 
		ymax = 3, yspace=10,nbins=50, overplot=True, track_agn = True, 
		track_subhalo_ptcls = False)
	print "zoom done"

	del(agn_cumul_ptcls)
	#del(subhalo_ptcls)# = None
	gc.collect()


#if __name__=="__main__":
def trace_outflows():
	#for snap in unique_snaps[57:]:
	for snapnum in [28, 41, 54,57]:#keysnaps: 
		#snapnum = unique_snaps.index(snap)
		snap = unique_snaps[snapnum]
		snap_halo = steps[snapnum].halos[halonum]
		stime = snap_halo.timestep.time_gyr
		print "t = %0.2f Gyr" % stime
		active_snap = pynbody.load(snap)
		halo_ptcls = active_snap.halos(dosort=True).load_copy(halonum_pb)
		halo_ptcls.physical_units()
		print "halo loaded"
		
		pynbody.analysis.halo.center(halo_ptcls, mode='ssc')
		#1: track particles from single burst
		if plot_rain:
		    rain(active_snap, halo_ptcls, snapnum)

		# # #2: track growing wind 
		if plot_zoom:
			zoom(halo_ptcls, steps, snapnum, halonum, snap, cumul_indices)

		print "snap %d done" % snapnum
		del(active_snap)
		del(halo_ptcls)
		gc.collect()


def mixing():
	hists = np.array([])
	for snap in unique_snaps[54:]:
		snapnum = unique_snaps.index(snap)
		snap_halo = steps[snapnum].halos[halonum]
		stime = snap_halo.timestep.time_gyr
		active_snap = pynbody.load(snap)
		halo_ptcls = active_snap.halos(dosort=True).load_copy(halonum_pb)
		halo_ptcls.physical_units()
		print "halo loaded"
		core_ptcls = halo_ptcls.g[(np.in1d(halo_ptcls.g['iord'], core_ind))]
		del(halo_ptcls)
		gc.collect()
		r = np.linalg.norm(core_ptcls['pos'] - snap_halo['shrink_center'], axis=1)
		hist, bins =np.histogram(r, range=(0,100), bins=100)
		del(r)
		gc.collect()
		colors= hist/float(sum(hist))
		hists = np.append(hists, colors)
		np.save('mixing_histogram', hists)

def plot_mixing():
		colors = np.load('mixing_histogram.npy')
		colors = np.reshape(colors, [len(colors)/100, 100])
		colors = np.log10(colors[colors > 0]) #log PDF of position at given snap
		print bins[hist != 0], colors
		plt.scatter(bins[:-1][hist!=0], [stime]*len(colors), color=cm.magma(colors))
		print 'snap %s done!' % snapnum
		plt.ylim(10,14)
		plt.savefig('mixing.png')
