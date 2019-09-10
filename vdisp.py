import numpy as np
import matplotlib
# if __name__ == "__main__":
#     matplotlib.use("Agg")
import matplotlib.pylab as plt
import pynbody, tangos, glob, gc 

sim = tangos.get_simulation('h1.cosmo50')
steps = sim.timesteps
stepname = [step.relative_filename for step in steps]
steptime = [step.time_gyr for step in steps]
stepnum = -1
halonum = 0 #generalize these later
halo = steps[stepnum].halos[halonum]

basename = '/nobackupp2/mtremmel/Romulus/'
datadir = basename+'h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096'

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

times = np.array([7.67, 9.38, 11.00, 11.65, 11.69, 11.97, 12.15, 13.8]) #
    
def vdisp(profile):
    return np.sqrt(profile['vx_disp']**2 + profile['vy_disp']**2 + profile['vz_disp']**2)

def get_profile(snap, weight='mass', rmin=5, rmax=2e3, tmax=None, tmin=None):
    if tmin:
        if tmax:
            ptcls = snap.g[(pynbody.filt.HighPass('temp',tmin))&(pynbody.filt.LowPass('temp',tmax))]
    elif tmax:
        ptcls = snap.g[pynbody.filt.LowPass('temp',tmax)]
    else:
        ptcls = snap.g
    return pynbody.analysis.profile.Profile(ptcls, min=rmin, max=rmax, type='log', ndim=3,weight_by=weight)

def vdisp_by_temp(times, profile=False, suffix='', rmin=5, rmax=2e3, tmin=1e6, tmax=1e6):
    h1 = halo
    while h1.timestep.time_gyr > times.min():
            h1 = h1.previous
    h1 = h1.next
    for time in times:
        while h1.timestep.time_gyr < time:
                h1 = h1.next
        stepind = np.argmin(abs(steptime - time))
        stepfile = unique_snaps[stepind]
        snap = pynbody.load(stepfile)
        snap.physical_units()
        print('snap collected')

        if profile:
            snap = snap[pynbody.filt.Sphere(h1['max_radius'], h1['shrink_center'])]
            print('sphere cutout')
        else:
            snap = snap[pynbody.filt.Annulus(0.03*h1.calculate('radius(500)'), 0.05*h1.calculate('radius(500)'), h1['shrink_center'])]
            print('ring cutout')
        
        snap['pos'] -= h1['shrink_center']
        snap['vel'] -= h1['Vcom']
        print('centered')
        
        if profile:
            pg = get_profile(snap, rmin = rmin, rmax =rmax)
            pghot = get_profile(snap, set='hot', rmin = rmin, rmax =rmax)
            pgcold = get_profile(snap, set='cold', rmin = rmin, rmax =rmax)
            print('profiles made')
        
            vdall = vdisp(pg)
            vdhot = vdisp(pghot)
            vdcold = vdisp(pgcold)
            print ('vdisp computed')
            
            plt.clf()
            plt.plot(range(len(vdall)), vdall, c='k', label='All')
            plt.plot(range(len(vdall)), vdhot, c='r', label=r'T > 10$^6$K')
            plt.plot(range(len(vdall)), vdcold, c='b', label=r'T < 10$^6$K')

            plt.legend()
            x = 10**np.linspace(np.log10(rmin),np.log10(rmax),100)
            xmin = np.argmin(abs(rmin - x))
            xmax = np.argmin(abs(rmax - x))
            plt.xlim(xmin,xmax)
            plt.ylim(1,1000)
            plt.savefig('vdisp_%0.2fGyr%s.png' % (time, suffix))
            del(snap, pg, pghot, pgcold)

        else:
            gas = np.sqrt(snap.g['v2'])
            hotgas = np.sqrt(snap.g[pynbody.filt.HighPass('temp', 1e6)]['v2'])
            coldgas = np.sqrt(snap.g[pynbody.filt.LowPass('temp', 1e6)]['v2'])

            hist, bins = np.histogram(gas, range=(gas.min(), gas.max()), bins = 100)
            histhot, bins = np.histogram(hotgas, range=(gas.min(), gas.max()), bins = 100)
            histcold, bins = np.histogram(coldgas, range=(gas.min(), gas.max()), bins = 100)

            plt.clf()
            plt.plot(bins[:-1], hist, c='k', label='All')
            plt.plot(bins[:-1], histhot, c='r', label=r'T > 10$^6$K')
            plt.plot(bins[:-1], histcold, c='b', label=r'T < 10$^6$K')
            plt.vlines(np.mean(gas), hist.min(), hist.max(), color='k')
            plt.vlines(np.median(gas), hist.min(), hist.max(), color='k', linestyle='dotted')
            plt.vlines(np.mean(hotgas), hist.min(), hist.max(), color='r')
            plt.vlines(np.median(hotgas), hist.min(), hist.max(), color='r', linestyle='dotted')
            plt.vlines(np.mean(coldgas), hist.min(), hist.max(), color='b')
            plt.vlines(np.median(coldgas), hist.min(), hist.max(), color='b', linestyle='dotted')

            plt.yscale('log')
            plt.legend()
            plt.xlim(0,5000)
            plt.savefig('vhist_%0.2fGyr.png' % time)
            plt.xlim(0,2000)
            plt.savefig('vhist_%0.2fGyr_xzoom.png' % time)
            del(snap, gas, hotgas, coldgas)
        
        gc.collect()
        print( "plots made")

if __name__=="__main__":
    vdisp_by_temp(times)

