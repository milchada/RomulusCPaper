#trace_outflows.py

# "NOTE THAT COLOR LIMITS DEPEND ON WHETHER IT'S A ZOOM"

import numpy as np
import tangos, pynbody, matplotlib, glob, gc 
# if __name__=="__main__":
#         matplotlib.use('Agg') #makes sure matplotlib doesn't look for a display when plottin
from tangos.examples import mergers
import matplotlib.pylab as plt
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

basename = '/nobackupp2/mtremmel/Romulus/'
datadir = basename+'h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096'


def plot(h1, ncoreptcl = 1000, width = 1700, suffix='', ent=True, pres=True, rho=False, temp=False, mach=False):
        merger_snap = h1.timestep.filename
        merger_sim = pynbody.load(merger_snap)
        "Merger snap loaded"
        merger_sim.physical_units()
        sph = merger_sim[pynbody.filt.Sphere(h1['max_radius'], h1['shrink_center'])]
        print "Sphere loaded "

        t2 = h1.timestep.time_gyr
        print "t = %.2f Gyr" % t2
        
        h1ptcls = sph.g
        del(sph, merger_sim)
        gc.collect()

        pynbody.analysis.halo.center(h1ptcls)
        print "centered"
        widthkpc = str(width)+' kpc'
        
        h1ptcls = h1ptcls[pynbody.filt.BandPass('z',"-1 Mpc","1 Mpc")] 
        
        h1ptcls['vdisp'] = np.sqrt(h1ptcls.g['vx_disp']**2 + h1ptcls.g['vy_disp']**2 + h1ptcls.g['vz_disp']**2)
        h1ptcls['vdisp'].units = 'km/s'
        
        h1ptcls['kT'] = h1ptcls.g['temp']*pynbody.units.k.in_units('keV K**-1')
        h1ptcls['kT'].units = 'keV'

        h1ptcls['entropy'] = h1ptcls.g['kT'] * pow(h1ptcls['ne']*h1ptcls.g['rho'].in_units('m_p cm**-3'), -2./3)
        h1ptcls['entropy'].units = 'keV cm**2'
        print "entropy calculated"
        h1ptcls['mach'] = np.sqrt(h1ptcls.g['v2'])/h1ptcls.g['cs']
        print "mach # calculated"

        h1ptcls.rotate_x(90) #so now x-z plane instead of x-y
        print "Rotated about x axis 90deg"
# if y:
        if ent:
                ent = image(h1ptcls, width=widthkpc,qty='entropy', qtytitle=r'K',
                title='%0.2f Gyr' % t2 , cmap=cm.magma, vmin=1, vmax=1e3)
                plt.xlim(-width/2,width/2)
                plt.ylim(-width/2,width/2)
                plt.savefig('sloshing_Kproj_xz_'+str(t2)[:5]+'_Gyr'+suffix+'.png' )
                plt.clf()
                del(ent)
                gc.collect()

        if pres:
                pressure = image(h1ptcls, width=widthkpc,qty='p', qtytitle=r'P', units='keV cm**-3',
                        title='%0.2f Gyr' % t2 , cmap=cm.viridis, vmin=1e-4, vmax=1)
                plt.xlim(-width/2,width/2)
                plt.ylim(-width/2,width/2)
                plt.savefig('sloshing_Pproj_xz_'+str(t2)[:5]+'_Gyr'+suffix+'.png' )
                plt.clf()
                del(pressure)
                gc.collect()

        if rho:
                rho = image(h1ptcls, width=widthkpc,qty='rho', qtytitle=r'$\rho$',
                        title='%0.2f Gyr' % t2 , cmap=cm.plasma, units='Msol kpc^-3', vmin=1e4, vmax=2e7)
                plt.xlim(-width/2,width/2)
                plt.ylim(-width/2,width/2)
                plt.savefig('sloshing_rhoproj_xz_'+str(t2)[:5]+'_Gyr'+suffix+'.png' )
                plt.clf()
                del(rho)
                gc.collect()

        if temp:
                temp = image(h1ptcls, width=widthkpc,qty='kT', qtytitle=r'T',
                        title='%0.2f Gyr' % t2 , cmap=cm.RdBu_r, vmin=0.1, vmax=5)
                plt.xlim(-width/2,width/2)
                plt.ylim(-width/2,width/2)
                plt.savefig('sloshing_Tproj_xz_'+str(t2)[:5]+'_Gyr'+suffix+'.png' )
                plt.clf()
                del(temp)
                gc.collect()

                mach = image(h1ptcls, width=widthkpc,qty='mach', qtytitle=r'M',
                        title='%0.2f Gyr' % t2 , cmap=cm.afmhot, vmin=0.01, vmax=10)
                plt.xlim(-width/2,width/2)
                plt.ylim(-width/2,width/2)
                plt.savefig('sloshing_Mach_xz_'+str(t2)[:5]+'_Gyr'+suffix+'.png' )
                plt.clf()
                del(mach, h1ptcls)
                gc.collect()
        print "xz finished"

def offset(merger=5, width = 1700, suffix=''):

        h1, h2 = mtree[2][merger]

        while h1.timestep.time_gyr > 9.76:
                h1 = h1.next

        plot(h1, width=width, suffix=suffix,ent=False, pres = False, rho=True, temp=True, mach=True)
        print "9.77 done"

        while h1.timestep.time_gyr < 10.99:
                h1 = h1.next
                
        plot(h1, width=width,ent=False, pres = False,  rho=True, temp=True, mach=True)
        print "11.00 Gyr done"

import glob
