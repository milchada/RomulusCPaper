import pynbody as pb 
import pynbody.plot.sph as sph
import matplotlib.pylab as plt

datadir = "/project/fas/nagai/lm643/h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096"
sim = pb.load(datadir)

h = sim.halos(dosort=True)

def plots(property='rho',firsthalo=0,nhalos=1):
	for halonum in xrange(firsthalo,nhalos):
		halo = h.load_copy(halonum)
		print "Halo %i loaded" % halonum
		halo.physical_units()
		halomass = sum(halo['mass'])
		halo.g['rhosq'] = halo.g['rho']*halo.g['rho']
		print "Density squared computed"
		#plot density profiles
		#bin by position
		pb.analysis.angmom.faceon(halo) #orient faceon wrt angular momentum
		print "Orientation complete"
			#skipping centering leads to IndexError in next steps
		volproj = sph.image(halo.g, qty=property, units='g cm^-3', av_z=True, width=100, cmap='RdYlBu',
							title = r"%e $M_\odot$" % halomass,
							filename='%i_volproj.png' % halonum) #slice
		print "Volume weighted projection complete"
		rhoproj = sph.image(halo.g, qty=property, units='g cm^-3', av_z='rhosq', width=100, cmap='RdYlBu',
							title = r"%e $M_\odot$" % halomass,
							filename='%i_emsnproj.png' % halonum) #integrated along z by default 
		print "Emission weighted projection complete" 
		# avg_sl, min_sl, max_sl, bins = sph.image_radial_profile(slice)
		# avg_pr, min_pr, max_pr, bins = sph.image_radial_profile(proj)

