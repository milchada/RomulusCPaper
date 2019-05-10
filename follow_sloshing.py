#follow sloshing
import numpy as np
import tangos
import pynbody
from tangos.examples import mergers
import matplotlib.pylab as plt
import gc
from pynbody.plot import image
from matplotlib import cm, colors

#load halo database
sim = tangos.get_simulation('h1.cosmo50')
steps = sim.timesteps
stepname = [step.relative_filename for step in steps]
steptime = [step.time_gyr for step in steps]
stepnum = -1
halonum = 0 #generalize these later
halo = steps[stepnum].halos[halonum]

mtree = mergers.get_mergers_of_major_progenitor(halo)

#major_mergers = [3,5] have mass ratios (and occur at) of 8.8 (11.69 Gyr) and 7.8 (11.32 Gyr), respectively

basename = '/nobackupp2/mtremmel/Romulus/'
datadir = basename+'h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096'

import glob
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

def plot(h1ptcls, h2ptcls, step):
	t2 = steptime[merger_ind + step]
	h1ptcls.physical_units()		
	h1pos = steps[merger_ind + step].halos[0]['shrink_center']
	
	h2ptcls.physical_units()
	phi = h2ptcls['phi']
	print "potentials gathered"
	sort = np.argsort(phi)
	print "most bound particles collected"
	del(phi)
	gc.collect()

	h1ptcls.g['entropy'] = h1ptcls.g['temp']*pynbody.units.k.in_units('keV K**-1') * pow(h1ptcls.g['rho'].in_units('m_p cm**-3'), -2./3)
	h1ptcls.g['entropy'].units = 'keV cm**2'
	print "entropy calculated"
	
	pynbody.analysis.halo.center(h1ptcls.g, mode='ssc')

	print "box centred on halo 1"
	ent = image(h1ptcls.g, width='500 kpc',qty='entropy', qtytitle=r'K', units='keV cm**2',  
                title='%0.2f Gyr' % t2, cmap=cm.magma, vmin=1, vmax=1e3)
	h2pos = h2ptcls['pos'][sort][:ncoreptcl]
	plt.scatter(h2pos[:,0],h2pos[:,1],alpha=0.15, c='w')#vel_colors, lw=0)
	plt.xlim(-250,250)
	plt.ylim(-250,250)
	plt.savefig('sloshing_proj_xy_%d_Gyr' % t2)
	print 'xy plane finished'

	h1ptcls.rotate_x(90) #so now x-z plane instead of x-y
	h2ptcls.rotate_x(90)
	print "Rotated about x axis 90º"
	ent = image(h1ptcls.g, width='1700 kpc',qty='entropy', qtytitle=r'K', units='keV cm**2',  
                title='%0.2f Gyr' % t2, cmap=cm.magma, vmin=1, vmax=1e3)
	h2pos = h2ptcls['pos'][sort][:ncoreptcl] 
	plt.scatter(h2pos[:,0],h2pos[:,1],alpha=0.15, c='w')#vel_colors, lw=0)
	plt.xlim(-850,850)
	plt.ylim(-850,850)
	plt.savefig('sloshing_proj_xz_%d_Gyr' % t2)	

	h1ptcls.rotate_y(90) #so now y-z plane instead of x-z
	h2ptcls.rotate_y(90)
	print "Rotated about x axis 90º"
	ent = image(h1ptcls.g, width='1700 kpc',qty='entropy', qtytitle=r'K', units='keV cm**2',  
                title='%0.2f Gyr' % t2, cmap=cm.magma, vmin=1, vmax=1e3)
	h2pos = h2ptcls['pos'][sort][:ncoreptcl] 
	plt.scatter(h2pos[:,0],h2pos[:,1],alpha=0.15, c='w')#vel_colors, lw=0)
	plt.xlim(-850,850)
	plt.ylim(-850,850)
	plt.savefig('sloshing_proj_yz_%d_Gyr' % t2)	

	
def offset(merger=5, ncoreptcl = 1000, endstep=4):
	h1, h2 = mtree[2][merger]
	
	merger_snap = str(basename+h2.path.split('/halo')[0])
	merger_sim = pynbody.load(merger_snap)
	"Merger snap loaded"
	h1ptcls = h1.load()
	print "Halo 1 loaded"
	h2ptcls = h2.load()
	print "Halo 2 loaded"
	h = merger_sim.halos()
	print "Halo cat loaded "
	merger_ind = steptime.index(h1.timestep.time_gyr)
	
	plot(h1ptcls, h2ptcls, 0)

	for step in xrange(1, endstep):
		plt.clf()
		current_snap = pynbody.load(unique_snaps[merger_ind+step])
		t2 = current_snap.properties['time'].in_units('Gyr') #tested - I'm using the right snapshots 
		print "Snap %d loaded" % step
		b = pynbody.bridge.OrderBridge(merger_sim, current_snap)
		print "Bridge made"
		
		h1ptcls = current_snap.halos(dosort=True).load_copy(1)
		del(current_snap)
		gc.collect()
		h2ptcls = b(h[h2.halo_number])
		print("particles collected")
		plot(h1ptcls, h2ptcls, step)
		del(h1ptcls, h2ptcls)
		gc.collect()
		print "Step %d done" % step

#omg i bet this is a projection issue