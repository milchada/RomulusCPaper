import glob
import numpy as np 
# import matplotlib
# if __name__=="__main__":
# 	matplotlib.use('Agg') #makes sure matplotlib doesn't look for a display when plotting
import matplotlib.pylab as plt
from matplotlib.pylab import cm
import gc
from scipy.stats import binned_statistic
import tangos as db
import pynbody 
import matplotlib

#load particles
basename = '/nobackupp2/mtremmel/Romulus/'
datadir = basename+'h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096'
# basename= '/gpfs/project/fas/nagai/etl28/Romulus/'
# datadir = basename+'RomulusC_Snapshots/h1.cosmo50PLK.1536gst1bwK1BH.004096'
sim = db.get_simulation('h1.cosmo50')
steps = sim.timesteps
halonum = 0
halo = steps[-1].halos[halonum]
steptime = halo.reverse_property_cascade('t()')[0]

keytimes = [7,8,10,11,11.3,11.65, 11.69, 11.97, 12.15] #Gyr
snapnums = [np.argmin(abs(np.array(steptime)-keytime)) for keytime in keytimes]

snapnum = 3
rcore = 10
radius = str(rcore)+' kpc'
nbins = 25

allgas = True

def filter(snap,halo_center=[0,0,0], radius='1 kpc'):
	snap.physical_units()
	spfilter = pynbody.filt.Sphere(radius, halo_center) & (pynbody.filt.FamilyFilter(pynbody.family.gas))# | pynbody.filt.FamilyFilter(pynbody.family.star))
	sphere = snap[spfilter]
	return sphere.g['iord']

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

def entropy(ptcls, allgas=False):
	T_kev = ptcls['temp'].in_units('K')*pynbody.units.k.in_units('keV K**-1')
	if allgas:
		n_cc = ptcls['rho'].in_units('m_p cm**-3')/ptcls['mu']
	else:
		n_cc = ptcls['ne']*ptcls.g['rho'].in_units('m_p cm**-3')
	entropy = T_kev*pow(n_cc,-2./3)
	entropy.units = 'keV cm^2'
	return entropy #*weight

# merger_sim = pynbody.load('/nobackupp2/mtremmel/Romulus/h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.003360')
def select_core_particles(snapnum, radius='30 kpc', halonum=1):
	sim = unique_snaps[snapnum]
	halo_ptcls = pynbody.load(sim).halos(dosort = True).load_copy(halonum)
	halo_ptcls.physical_units()
	pynbody.analysis.halo.center(halo_ptcls, mode='ssc')
	core_ind = filter(halo_ptcls, radius=radius)
	del(halo_ptcls)# = None #clear cache
	gc.collect()
	return core_ind

def match_gas_ind(original_core_ind, step_core_ptcls):
	index = np.argsort(original_core_ind)
	sorted_x = original_core_ind[index]
	sorted_index = np.searchsorted(sorted_x, step_core_ptcls['iord'])
	yindex = np.take(index, sorted_index, mode="clip")
	del(index, sorted_index, sorted_x)
	gc.collect()
	mask = original_core_ind[yindex] != step_core_ptcls['iord']
	xind = np.ma.array(yindex, mask=mask).data
	del(mask, yindex)
	gc.collect()
	return xind 
	
def collect_ptcls(filename):
	core_ind = select_core_particles(snapnums[snapnum - 1], radius)
	Ks = np.empty([len(core_ind), 4])
	halo_ptcls = pynbody.load(unique_snaps[snapnums[snapnum - 1]]).halos(dosort=True).load_copy(1)
	core_gas = halo_ptcls.g[(np.in1d(halo_ptcls.g['iord'], core_ind))]
	del(halo_ptcls)
	gc.collect()
	Ki = entropy(core_gas, allgas=allgas)
	del(core_gas)
	gc.collect()
	Ks[:,0] = core_ind
	Ks[:,1] = Ki
	del(core_ind, Ki)
	gc.collect()

	halo_ptcls = pynbody.load(unique_snaps[snapnums[snapnum]]).halos(dosort=True).load_copy(1)
	core_gas = halo_ptcls.g[(np.in1d(halo_ptcls.g['iord'], Ks[:,0]))]
	del(halo_ptcls)
	gc.collect()
	K = entropy(core_gas,allgas=allgas)
	mass = core_gas['mass'].in_units('Msol')
	gas_ind = match_gas_ind(Ks[:,0], core_gas) #Ks[:,0]
	Ks[gas_ind, 2] = K
	Ks[gas_ind, 3] = mass
	del(core_gas, K, gas_ind)
	gc.collect()
	for ptcl in xrange(len(Ks)):
		if Ks[ptcl, 2] < 1e-5: 	#get rid of empty particles - they cause the scatter blow up and mean to fall
			Ks[ptcl,2] = np.nan
	np.save(filename, Ks)

def plot(Ks, tmin, tmax, ret_cbar=False, norm=None, m=None):
	binned = binned_statistic(np.log10(Ks[:,1]), np.log10(Ks[:,2]), bins = nbins, range = (-2, 2))
	x = 10**binned.bin_edges[:-1]
	median = []
	min = []
	max = []
	mcells = []
	for i in xrange(nbins):
		median.append(np.nanpercentile(Ks[:,2][binned.binnumber == i], 50))
		min.append(np.nanpercentile(Ks[:,2][binned.binnumber == i], 25))
		max.append(np.nanpercentile(Ks[:,2][binned.binnumber == i], 75))
		mcells.append(sum(Ks[:,3][binned.binnumber == i]))
	median = np.array(median)
	min = np.array(min)
	max = np.array(max)
	mcells = np.array(mcells)
	plt.clf()
	# plt.plot(x, x, c='g') #ax[row,col]
	plt.plot(x, median, c='w')
	# plt.fill_between(x, min, max, color=colors,alpha=0.5)
	fig, ax = plt.subplots()
	fig1, ax1 = plt.subplots()
	rand = np.random.random_integers(low = mcells.min(), high=mcells.max(), size=[50,50])
	if ret_cbar:
		norm = matplotlib.colors.LogNorm(vmin = mcells.min(),vmax = mcells.max())
		m = cm.ScalarMappable(norm = norm, cmap = cm.RdBu)
	f1 = ax.imshow(rand, cmap = cm.RdBu, norm=norm)
	ax1.bar(x, bottom=min, height=(max - min), width=10**(binned.bin_edges[1:])-x ,color=m.to_rgba(mcells))
	ax1.plot(x, x, c='g')
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax1.xaxis.set_tick_params(labelsize='x-large')
	ax1.yaxis.set_tick_params(labelsize='x-large')
	ax1.set_ylim(1e-3, 1e2)
	ax1.set_xlim(1e-2, 1e2)
	ax1.set_xlabel(r'K (t$_i$)', fontsize = 'x-large' )
	ax1.set_ylabel(r'K (t$_f$)', fontsize = 'x-large' )
	fig.colorbar(f1, ax=ax1)
	plt.tight_layout()
	print( "plotting complete")
	if allgas:
		suffix = '_allgas'
	else:
		suffix = ''
	plt.savefig('mixing_%0.2fGyr_%dkpc%s.png' % (tmin, rcore, suffix))
	if ret_cbar:
		return norm, m

	# #consider particles within sphere, not only in halo
	# #what fraction am i losing? maybe track mass as well
	# #send around updated draft 
def plotall():
	import glob
	files = glob.glob('*Kbins*.npy')
	files.sort()

	norm, m = plot(np.load(files[0]), 10, 11, ret_cbar = True)
	plot(np.load(files[1]), 11, 11.3, norm=norm, m=m)
	plot(np.load(files[3]), 11.65, 11.69, norm=norm, m=m)
	plot(np.load(files[5]), 11.97, 12.15, norm=norm, m=m)

	#plot the 200kpc ones too
	
	# norm, m = plot(np.load(files[4]), 10, 11, ret_cbar = True)
	# plot(np.load(files[5]), 11, 11.3, norm=norm, m=m)
	# plot(np.load(files[6]), 11.65, 11.69, norm=norm, m=m)
	# plot(np.load(files[7]), 11.97, 12.15, norm=norm, m=m)

