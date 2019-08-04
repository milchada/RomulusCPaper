#trace_outflows.py

import numpy as np
import tangos
import pynbody
import matplotlib
if __name__=="__main__":
        matplotlib.use('Agg') #makes sure matplotlib doesn't look for a display when plottin
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

def plot(h1ptcls, merger_ind, step, ncoreptcl = 1000, x=True, y=True, z=True, width = 1700, subhalo=False, h2ptcls=None, suffix=''):
        widthkpc = str(width)+' kpc'
        #for test only 
        ncoreptcl = 1000
        t2 = steptime[merger_ind + step]
        print "t = %.2f Gyr" % t2
        
        h1ptcls = h1ptcls[pynbody.filt.BandPass('z',"-1 Mpc","1 Mpc")] 
        if subhalo:
                h2pos = h2ptcls['pos']
                print "positions gathered"
                sort = np.argsort(np.linalg.norm(h2pos, axis=1))
                h2pos = h2ptcls['pos'][sort][:ncoreptcl]
        print "particles closest to halo 1 collected"

        h1ptcls['kT'] = h1ptcls.g['temp']*pynbody.units.k.in_units('keV K**-1')
        h1ptcls['kT'].units = 'keV'

        h1ptcls['entropy'] = h1ptcls.g['kT'] * pow(h1ptcls['ne']*h1ptcls.g['rho'].in_units('m_p cm**-3'), -2./3)
        h1ptcls['entropy'].units = 'keV cm**2'
        print "entropy calculated"

        if x:
                ent = image(h1ptcls, width=widthkpc,qty='entropy', qtytitle=r'K',
                                title='%0.2f Gyr' % t2, cmap=cm.magma, vmin=1, vmax=1e3)
                if subhalo:
                        plt.scatter(h2pos[:,0],h2pos[:,1],alpha=0.15, c='w')
                plt.xlim(-width/2,width/2)
                plt.ylim(-width/2,width/2)
                plt.savefig('sloshing_Kproj_xz_%0.2f_Gyr'+suffix+'.png' % t2)
                plt.clf()
                del(ent)
                gc.collect()

                pressure = image(h1ptcls, width=widthkpc,qty='p', qtytitle=r'P', units='keV cm**-3',
                        title='%0.2f Gyr' % t2, cmap=cm.viridis, vmin=1e-4, vmax=1)
                if subhalo:
                        plt.scatter(h2pos[:,0],h2pos[:,1],alpha=0.15, c='w')
                plt.xlim(-width/2,width/2)
                plt.ylim(-width/2,width/2)
                plt.savefig('sloshing_Pproj_xz_%0.2f_Gyr'+suffix+'.png' % t2)
                plt.clf()
                del(pressure)
                gc.collect()

                rho = image(h1ptcls, width=widthkpc,qty='rho', qtytitle=r'$\rho$',
                        title='%0.2f Gyr' % t2, cmap=cm.plasma, units='Msol kpc^-3', vmin=5e3, vmax=5e7)
                if subhalo:
                        plt.scatter(h2pos[:,0],h2pos[:,1],alpha=0.15, c='w')
                plt.xlim(-width/2,width/2)
                plt.ylim(-width/2,width/2)
                plt.savefig('sloshing_rhoproj_xz_%0.2f_Gyr'+suffix+'.png' % t2)
                plt.clf()
                del(rho)
                gc.collect()

                temp = image(h1ptcls, width=widthkpc,qty='kT', qtytitle=r'T',
                        title='%0.2f Gyr' % t2, cmap=cm.RdBu_r, vmin=0.1, vmax=10)
                if subhalo:
                        plt.scatter(h2pos[:,0],h2pos[:,1],alpha=0.15, c='w')
                plt.xlim(-width/2,width/2)
                plt.ylim(-width/2,width/2)
                plt.savefig('sloshing_Tproj_xz_%0.2f_Gyr'+suffix+'.png' % t2)
                plt.clf()
                del(temp)
                gc.collect()


        if y or z:
                h1ptcls.rotate_x(90) #so now x-z plane instead of x-y
                print "Rotated about x axis 90deg"
        if y:
                filename = 'sloshing_proj_xz_'+str(t2)+'_Gyr'+suffix+'.png'
                if not glob.glob(filename):
                        ent = image(h1ptcls, width=widthkpc,qty='entropy', qtytitle=r'K',
                                title='%0.2f Gyr' % t2, cmap=cm.magma, vmin=1, vmax=1e3)
                        if subhalo:
                                plt.scatter(h2pos[:,0],h2pos[:,2],alpha=0.15, c='w')
                        plt.xlim(-width/2,width/2)
                        plt.ylim(-width/2,width/2)
                        plt.savefig('sloshing_Kproj_xz_%0.2f_Gyr'+suffix+'.png' % t2)
                        plt.clf()
                        del(ent)
                        gc.collect()

                        pressure = image(h1ptcls, width=widthkpc,qty='p', qtytitle=r'P', units='keV cm**-3',
                                title='%0.2f Gyr' % t2, cmap=cm.viridis, vmin=1e-4, vmax=1)
                        if subhalo:
                                plt.scatter(h2pos[:,0],h2pos[:,2],alpha=0.15, c='w')
                        plt.xlim(-width/2,width/2)
                        plt.ylim(-width/2,width/2)
                        plt.savefig('sloshing_Pproj_xz_%0.2f_Gyr'+suffix+'.png' % t2)
                        plt.clf()
                        del(pressure)
                        gc.collect()

                        rho = image(h1ptcls, width=widthkpc,qty='rho', qtytitle=r'$\rho$',
                                title='%0.2f Gyr' % t2, cmap=cm.plasma, units='Msol kpc^-3', vmin=5e3, vmax=5e7)
                        if subhalo:
                                plt.scatter(h2pos[:,0],h2pos[:,2],alpha=0.15, c='w')
                        plt.xlim(-width/2,width/2)
                        plt.ylim(-width/2,width/2)
                        plt.savefig('sloshing_rhoproj_xz_%0.2f_Gyr'+suffix+'.png' % t2)
                        plt.clf()
                        del(rho)
                        gc.collect()

                        temp = image(h1ptcls, width=widthkpc,qty='kT', qtytitle=r'T',
                                title='%0.2f Gyr' % t2, cmap=cm.RdBu_r, vmin=0.1, vmax=10)
                        if subhalo:
                                plt.scatter(h2pos[:,0],h2pos[:,2],alpha=0.15, c='w')
                        plt.xlim(-width/2,width/2)
                        plt.ylim(-width/2,width/2)
                        plt.savefig('sloshing_Tproj_xz_%0.2f_Gyr'+suffix+'.png' % t2)
                        plt.clf()
                        del(temp)
                        gc.collect()

                        print "xz finished"

        if z:
                filename = 'sloshing_proj_xz_'+str(t2)+'_Gyr'+suffix+'.png'
                if not glob.glob(filename):
                        h1ptcls.rotate_y(90) #so now y-z plane instead of x-z
                        print "Rotated about y axis 90deg"
                        
                        ent = image(h1ptcls, width=widthkpc,qty='entropy', qtytitle=r'K',
                                title='%0.2f Gyr' % t2, cmap=cm.magma, vmin=1, vmax=1e3)
                        if subhalo:
                                plt.scatter(h2pos[:,1],h2pos[:,2],alpha=0.15, c='w')
                        plt.xlim(-width/2,width/2)
                        plt.ylim(-width/2,width/2)
                        plt.savefig('sloshing_Kproj_yz_%0.2f_Gyr'+suffix+'.png' % t2)
                        plt.clf()
                        del(ent)
                        gc.collect()

                        pressure = image(h1ptcls, width=widthkpc,qty='p', qtytitle=r'P', units='keV cm**-3',
                                title='%0.2f Gyr' % t2, cmap=cm.viridis, vmin=1e-4, vmax=1)
                        if subhalo:
                                plt.scatter(h2pos[:,1],h2pos[:,2],alpha=0.15, c='w')
                        plt.xlim(-width/2,width/2)
                        plt.ylim(-width/2,width/2)
                        plt.savefig('sloshing_Pproj_yz_%0.2f_Gyr'+suffix+'.png' % t2)
                        plt.clf()
                        del(pressure)
                        gc.collect()

                        rho = image(h1ptcls, width=widthkpc,qty='rho', qtytitle=r'$\rho$',
                                title='%0.2f Gyr' % t2, cmap=cm.plasma, units='Msol kpc^-3', vmin=5e3, vmax=5e7)
                        if subhalo:
                                plt.scatter(h2pos[:,1],h2pos[:,2],alpha=0.15, c='w')
                        plt.xlim(-width/2,width/2)
                        plt.ylim(-width/2,width/2)
                        plt.savefig('sloshing_rhoproj_yz_%0.2f_Gyr'+suffix+'.png' % t2)
                        plt.clf()
                        del(rho)
                        gc.collect()

                        temp = image(h1ptcls, width=widthkpc,qty='kT', qtytitle=r'T',
                                title='%0.2f Gyr' % t2, cmap=cm.RdBu_r, vmin=0.1, vmax=10)
                        if subhalo:
                                plt.scatter(h2pos[:,1],h2pos[:,2],alpha=0.15, c='w')
                        plt.xlim(-width/2,width/2)
                        plt.ylim(-width/2,width/2)
                        plt.savefig('sloshing_Tproj_yz_%0.2f_Gyr'+suffix+'.png' % t2)
                        plt.clf()
                        del(temp)
                        gc.collect()

                        print "yz finished"

def offset(merger=5, tmin = 11.64, startstep = 0, endstep=4, width = 1700, subhalo=False, suffix=''):

        h1, h2 = mtree[2][merger]
        while h1.timestep.time_gyr < tmin:
                h1 = h1.next
                h2 = h2.next
        merger_snap = str(basename+h2.path.split('/halo')[0])
        merger_sim = pynbody.load(merger_snap)
        "Merger snap loaded"
        h = merger_sim.halos()
        print "Halo cat loaded "
        merger_ind = steptime.index(h1.timestep.time_gyr)

        if not startstep:
                h1ptcls = h1.load().g
                h1ptcls.physical_units()
                # h2ptcls = h2.load().g
                # h2ptcls.physical_units()

                h1ptcls['pos'] -= h1['shrink_center']
                # h2ptcls['pos'] -= h1['shrink_center']
                
                plot(h1ptcls, merger_ind, 0, x=False, width=width, subhalo=subhalo, suffix=suffix)

        i= 0
        while i < startstep:
                h1 = h1.next
                i += 1

        for step in xrange(max(startstep,1), endstep):

                plt.clf()
                current_snap = pynbody.load(unique_snaps[merger_ind+step])
                print "Snap %d loaded" % step
                

                h1ptcls = h1.load().g
                h1ptcls.physical_units()
                h1ptcls['pos'] -= h1['shrink_center']

                if subhalo:
                        b = pynbody.bridge.OrderBridge(merger_sim, current_snap)
                	print "Bridge made"
                        h2ptcls = b(h[h2.halo_number])
                        print("particles collected")
                        h2ptcls.physical_units()
                        h2ptcls['pos'] -= h1['shrink_center']
                	print( "particles centered on halo 1")
                        plot(h1ptcls, merger_ind, step, x=False, width=width, subhalo=True, h2ptcls=h2ptcls)
                        del(current_snap, h1ptcls, h2ptcls, b)
                else:
                        plot(h1ptcls, merger_ind, step, x=False, width=width)
                        del(current_snap, h1ptcls)
                
                gc.collect()
                print "Step %d done" % step
                h1 = h1.next

if __name__=="__main__":
        offset(5, 11.64, startstep=4)
        offset(5, 11.64, startstep=0,width=200, suffix='zoom')
