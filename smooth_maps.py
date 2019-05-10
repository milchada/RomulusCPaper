import glob
# import pynbody
# from pynbody.plot import image
from astropy import units, cosmology, convolution
import matplotlib.pylab as plt
from matplotlib import cm
import matplotlib
import numpy as np
import yt
from yt.units import mp

#create emission-weighted bulk v_los, sigma_v maps projected in 2d

fov = {'xarm':3*60, 'athena':4*60, 'lynx':5*60} #arcsec; Athena is really 5' hexagon but i'm converting to equivalent square
pixel_side = {'xarm':6, 'athena': 62, 'lynx': 5*60} #pixel/side
psf = {'xarm':60, 'athena':5, 'lynx':.5} #arcsec, this is FWHM
spec_res = {'xarm':6, 'athena':2.5, 'lynx':3} #eV

lcdm = cosmology.default_cosmology.get()

###FOR ROMULUSC
basename = '/nobackupp2/mtremmel/Romulus/'
datadir = basename+'h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096'
ds = yt.load(datadir)
R200 = ds.quan(1032.7, "kpc")

z_perseus = 0.01790
h = pynbody.load(datadir).halos(dosort = True).load_copy(1)
print "halo loaded"
pynbody.analysis.angmom.faceon(h)
print "sim box oriented"
h.physical_units()
h.g['rhosq'] = h.g['rho']**2

def observation(instrument, redshift, filename, quantity = 'v_disp', qtytitle=r'$\sigma$', 
###END ROMULUSC

###FOR TNG
# ds = yt.load('/nobackupp8/jzuhone/tng/halo_3.hdf5')

# # Center of cluster
# C = ds.arr([127713.06189105, 120893.9511367 , 77686.08604591], "kpc")
# # R200 of cluster
# R200 = ds.quan(2035.5552580985202, "kpc")

# def _emission_measure(field,data):
# 	return data[(field,'rho')]**2

# ds.add_field(("gas","rhosq"), function=_emission_measure)

# def _emission_measure(field, data):
#     nenh = data["PartType0","Density"]*data["PartType0",'particle_mass']
#     nenh /= mp*mp
#     nenh.convert_to_units("cm**-3")
#     X_H = 0.76
#     nenh *= X_H * data["PartType0", 'ElectronAbundance']
#     nenh *= X_H * (1.-data["PartType0", 'NeutralHydrogenAbundance'])
#     return nenh
# ds.add_field(("PartType0", 'emission_measure'),
#              function=_emission_measure,
#              particle_type=True,
#              units="cm**-3", force_override=True)

# def _vz_squared(field, data):
#     return data["gas","velocity_z"]*data["gas","velocity_z"]
# ds.add_field(("gas","velocity_z_squared"), _vz_squared, units="cm**2/s**2", force_override=True)

# This defines a box around the cluster with width 4.0*R200 just to get everything
# le = C-2.0*R200
# re = C+2.0*R200
# reg = ds.box(le, re)

# # This projects along the line of sight
# prj = reg.integrate(("gas","velocity_z"), weight=("PartType0","emission_measure"), axis="z")
# ###END TNG

# def observation(instrument, redshift, filename, qtytitle=r'$\sigma$', 
	width=5000, vmin=1, vmax=1000, suffix='', step=150):
	Mpc_rad = lcdm.angular_diameter_distance(redshift)
	kpc_arcsec = Mpc_rad.to('kpc')/units.radian.in_units('arcsec')
	width_arcsec = width/kpc_arcsec #arcsec/side 
	pixel_res = fov[instrument]/pixel_side[instrument] #arcsec/pix for instrument
	resolution = (width_arcsec/pixel_res).value

	# nx = int(resolution)
	# frb = prj.to_frb((5.0, "Mpc"), nx, center=C)

	# Compute sigma2
	# sigma2 = frb["gas","velocity_z_squared"].to_value("km**2/s**2")-frb["gas","velocity_z"].to_value("km/s")**2
	# # Compute emission measure
	# EM = frb["PartType0","emission_measure"].d
	# #block/bin image same as instrument
	vdisp = image(h.g, resolution = resolution, qty=quantity, width=str(width)+' kpc', qtytitle=qtytitle, 
		av_z='rhosq', title='%0.2f' % redshift, cmap=cm.magma, vmin=1, vmax=max(250,vmax), log=False, noplot=True)
	print("image binned")

	EM = image(h.g, resolution = resolution, qty='rhosq', width=str(width)+' kpc', qtytitle=qtytitle, 
		av_z='vol', title='%0.2f' % redshift, cmap=cm.magma, vmin=1, vmax=max(250,vmax), log=False, noplot=True)

	fwhm_pix = psf[instrument]/pixel_res
	std_pix = fwhm_pix/2.355
	kernel = convolution.Gaussian2DKernel(stddev = std_pix)
	conv = convolution.convolve(np.sqrt(vdisp)*rhosq, kernel)/convolution.convolve(EM, kernel)
	print("image convolved with PSF")
	xticks = np.linspace(-width/2., width/2., int(resolution))/kpc_arcsec.value
	fig, ax = plt.subplots()
	a1=ax.imshow(conv.T, cmap = cm.magma, norm=matplotlib.colors.Normalize(vmin,vmax))
	c = (sigma2.shape[0]/2, sigma2.shape[0]/2)
	r500_arcsec = 1400/kpc_arcsec.value #change 1400 to r500 in kpc
	rad =r500_arcsec/pixel_res
	circ = matplotlib.patches.Circle(c, rad,edgecolor='w',linestyle='dotted', linewidth=3,fill=False)
	ax.add_patch(circ)
	rad /= 0.7 #r500 --> r200
	circ = matplotlib.patches.Circle(c, rad,edgecolor='w',linestyle='dotted', linewidth=3,fill=False)
	ax.add_patch(circ)
	plt.xticks(np.arange(int(resolution)/step)*step, ['%d' % tick for tick in xticks[::step]])
	plt.yticks(np.arange(int(resolution)/step)*step, ['%d' % tick for tick in xticks[::step]])
	plt.xlabel('R (")')
	plt.ylabel('R (")')
	plt.colorbar(a1,ax=ax)
	plt.savefig(filename+'_'+instrument+suffix)

if __name__=="__main__":
 	
 	observation('lynx', 0.3, 'snapz0_z03',vmax=200,suffix='zoom')
 	observation('athena', 0.3, 'snapz0_z03',vmax=75,suffix='zoom')
