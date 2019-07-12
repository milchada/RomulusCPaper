
import numpy as np 
from scipy.stats import binned_statistic
import tangos as db
import glob, gc, pynbody

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

keytimes = [9.2,9.7, 11.3, 11.65, 11.97, 12.20, 12.50,12.80, 13.25, 13.55] #Gyr
snapnums = [np.argmin(abs(np.array(steptime)-keytime)) for keytime in keytimes]
snapinds = [1,3,5,7,9] #
rcore = 30
radius = str(rcore)+' kpc'
nbins = 25

allgas = True
nosubs = True

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
		n_cc = ptcls['ne']*ptcls.g['rho'].in_units('m_p cm**-3')
	else:
		n_cc = ptcls.g['rho'].in_units('m_p cm**-3')/ptcls.g['mu']
	entropy = T_kev*pow(n_cc,-2./3)
	entropy.units = 'keV cm^2'
	return entropy #*weight

# merger_sim = pynbody.load('/nobackupp2/mtremmel/Romulus/h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.003360')
def select_core_particles(snapnum, radius='30 kpc', halonum=1):
	sim = unique_snaps[snapnum]
	halo_ptcls = pynbody.load(sim).halos(dosort = True).load_copy(halonum)
	if nosubs:
		halo_ptcls = halo_ptcls[halo_ptcls['amiga.grp'] == 1]
	halo_ptcls.physical_units()
	pynbody.analysis.halo.center(halo_ptcls, mode='ssc')
	core_ind = filter(halo_ptcls, radius=radius)
	core_gas = halo_ptcls.g[(np.in1d(halo_ptcls.g['iord'], core_ind))]
	del(halo_ptcls)# = None #clear cache
	gc.collect()
	return core_ind, core_gas

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
	
def collect_ptcls(snapnum,filename):
	core_ind, core_gas = select_core_particles(snapnums[snapnum-1], radius)
	core_gas.physical_units()
	print "halo particles loaded"
	Ks = np.empty([len(core_ind), 6]) #this array stores the properties of the particles that are in the core at time t1
	# Ks.dtype.names = ['ind','mass','Ki','Kf','ri','rf']
	Ks[:,0] = core_ind
	Ks[:,1] = core_gas['mass'].in_units('Msol')
	Ks[:,2] = entropy(core_gas, allgas=allgas)
	print "entropy computed"
	del(core_ind)
	gc.collect()
	Ks[:,4] = np.linalg.norm(core_gas['pos'], axis=1)
	print "dist computed"
	del(core_gas)
	gc.collect()

	#load the next snapshot

	halo_ptcls = pynbody.load(unique_snaps[snapnums[snapnum]]).halos(dosort=True).load_copy(1)
	if nosubs:
		halo_ptcls = halo_ptcls[halo_ptcls['amiga.grp'] == 1]
	halo_ptcls.physical_units()
	print "next snap loaded"
	pynbody.analysis.halo.center(halo_ptcls, mode='ssc')
	print "snap centred"
	core_gas = halo_ptcls.g[(np.in1d(halo_ptcls.g['iord'], Ks[:,0]))]
	del(halo_ptcls)
	gc.collect()
	
	gas_ind = match_gas_ind(Ks[:,0], core_gas) #Ks[:,0]
	Ks[gas_ind, 3] = entropy(core_gas, allgas=allgas)
	print "entropy computed"
	Ks[gas_ind, 5] = np.linalg.norm(core_gas['pos'], axis=1)
	print "dist computed"
	del(core_gas, gas_ind)
	gc.collect()
	for ptcl in xrange(len(Ks)):
		if Ks[ptcl, 2] < 1e-5:  #get rid of empty particles - they make the scatter blow up and mean to fall
			Ks[ptcl] = np.nan
	print "done!"
	np.save(filename, Ks)

if __name__=="__main__":
	from multiprocessing import Process
	proc = []
	for snapind in snapinds:
		filename = '%0.2f-%0.2f_Ks_histogram_%d_kpc' % (keytimes[snapind-1], keytimes[snapind], rcore)
		print snapind, filename
		p = Process(target = collect_ptcls, args=(snapind, filename))
		p.start()
		proc.append(p)
	for p in proc:
		p.join()