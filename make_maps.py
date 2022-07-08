import glob
import numpy as np 
import matplotlib
if __name__=="__main__":
	matplotlib.use('Agg') #makes sure matplotlib doesn't look for a display when plotting

import matplotlib.pylab as plt
import pynbody 
from pynbody.plot import image
from matplotlib.pylab import cm
import gc

#load particles
basename = '/nobackupp2/mtremmel/Romulus/'
datadir = basename+'h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096'

#load halo database
halonum_pb = 1

def entropy(ptcls):
	T_kev = ptcls['temp'].in_units('K')*pynbody.units.k.in_units('keV K**-1')
	n_cc = ptcls['ne']*ptcls.g['rho'].in_units('m_p cm**-3')
	entropy = T_kev*pow(n_cc,-2./3)
	# entropy.units = 'keV cm^2'
	return entropy #*weight

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
#bh_activity()

def maps(maptype, select, stime, snapnum, width='100 kpc', retimg = True):
	if maptype == 'rhoproj':
		img = image(select.g, width=width, qtytitle=r'$\rho$',units='Msol kpc^-2', 
		title='%0.2f' % stime, cmap=cm.viridis, vmin=1e5, vmax=1e9,filename='merger_density_proj_%d' % snapnum)

	if maptype == 'rho':
		img = image(select.g, width=width, qtytitle=r'$\rho$',units='Msol kpc^-3', 
		av_z='mass', title='%0.2f' % stime, cmap=cm.magma, vmin=5e3, vmax=5e7, filename='merger_density_%d' % snapnum)
	
	if maptype == 'temp':
		img = image(select.g, qty='temp',width=width, qtytitle='temp', units='K', log='False',
			av_z='mass', title='%0.2f' % stime, cmap=cm.RdBu_r, vmin = 1e6, vmax = 5e7, filename='merger_temp_%d' % snapnum)
	
	if maptype == 'entropy':
		img = image(select.g, qty='entropy',width=width, qtytitle='K', units='keV cm**2', log='True',
			av_z='mass', title='%0.2f' % stime, cmap=cm.magma, vmin = 1, vmax = 1e3, filename='merger_entropy_%d' % snapnum)
	
	if retimg:
		return img

keysnaps = [30,38,46,54,68]
colors = {30:0.05, 38:0.2, 46:0.4, 54:0.6, 68:0.8}
labels = {30:'Pre-quenching', 38: 'Quenching', 46:'Quenched', 54:'Pre-disruption', 68:'Disrupting'}


def make_map(snapnum, maptypes=['rho'],width='500 kpc'):
	snap = unique_snaps[snapnum]
	active_snap = pynbody.load(snap)
	# active_snap.physical_units()
	stime = active_snap.properties['time'].in_units('Gyr')
	print( "t = %0.2f Gyr" % stime)
	# halo_ptcls = active_snap.load_copy()#halos(dosort=True).load_copy(halonum_pb)
	# halo_ptcls.physical_units()
	print( "halo loaded")
	
	#3: zoom on halo, overplot mass-weighted average 
	# pynbody.analysis.halo.center(halo_ptcls.g, mode='com') #gas centroid
	# pynbody.analysis.halo.center(halo_ptcls.g, mode='ssc')
	# pynbody.analysis.angmom.faceon(halo_ptcls.g)
	# pynbody.analysis.halo.center(halo_ptcls.g, mode='ssc')
	halo_ptcls = active_snap
	slicethickness = 5
	select = halo_ptcls.g[pynbody.filt.BandPass('z',-1*slicethickness, slicethickness)]# & pynbody.filt.LowPass('temp',1e5)]
	select.physical_units()
	# select['mdotr'] = select['vr']*select['mass']/select['r']
	del(halo_ptcls)
	# mdotslice=image(select.g, qty='mdotr',width='100 kpc', qtytitle=r'$\dot{M}$', units='Msol yr**-1',
	# 	title='%0.2f' % stime, cmap=cm.coolwarm,vmin=-0.08, vmax= 0.08, log=False, filename='cold_mass_inflow_5kpc_slice_%d' % snapnum)

	# fig, ax = plt.subplot()
	# pynbody.plot.sph.velocity_image(select, width='100 kpc', cmap = "viridis", mode='stream', units='Msol kpc^-2',
	#                    density = 2.0, vector_resolution=100, vmin=1e4, vmax=5e8,subplot=ax,
	#                    show_cbar=False, vector_color='white')
	# fig.savefig('coldgas_v_%d.png' % snapnum)		
	# select.g['T_keV'] = select.g['temp'].in_units('K')*pynbody.units.k.in_units('keV K**-1')
	# select.g['T_keV'].units='keV'
	# # select['entropy'] = entropy(select)
	
	for maptype in maptypes:
		maps(maptype, select, stime, snapnum, width=width)

	del(select)
	gc.collect()
	
	print ("snap %d done" % snapnum)

if __name__=="__main__":
	for snapnum in keysnaps[1:]:
		make_map(snapnum, maptype='rho')