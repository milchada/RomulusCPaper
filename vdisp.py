import numpy as np
import pynbody, tangos, glob, gc 

sim = tangos.get_simulation('h1.cosmo50')
steps = sim.timesteps
stepnum = -1
halonum = 0 #generalize these later
halo = steps[stepnum].halos[halonum]

steps = [0, 8, 9, 10, 11, 12, 16, 17, 18, 25, 26, 27, 37, 38, 39]
ts = halo.calculate_for_progenitors('t()')[0]
times = np.array([ts[i] for i in steps])

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


def write(snap, i, c, vcom, sigma, sigmar, suffix, steps=steps, rmin=1, rmax=1e3):
    snap['pos'] -= c
    snap['vel'] -= vcom

    print('centered')

    pg = get_profile(snap.g, rmin = rmin, rmax =rmax)
    print('profiles made')

    sigma[steps[i]] = vdisp(pg)
    print ('vdisp computed')
    sigmar[steps[i]] = pg['vr_disp']
    np.save('sigma3d%s' % suffix, sigma)
    np.save('sigma_r%s' % suffix, sigmar)

    snap['pos'] += c
    snap['vel'] += vcom
    del(pg)
    gc.collect()

sigma1 = np.load('sigma3d_cdb_vcdb.npy')#np.zeros((length,100)) 
sigmar1 = np.load('sigma_r_cdb_vcdb.npy')#np.zeros((length,100)) 
sigma2 = np.load('sigma3d_cdb_vcom.npy')#np.zeros((length,100)) 
sigmar2 = np.load('sigma_r_cdb_vcom.npy')#np.zeros((length,100)) 
sigma3 = np.load('sigma3d_ccom_vcom.npy')#np.zeros((length,100)) 
sigmar3 = np.load('sigma_r_ccom_vcom.npy')#np.zeros((length,100)) 
sigma4 = np.load('sigma3d_ccom_vcdb.npy')#np.zeros((length,100)) 
sigmar4 = np.load('sigma_r_ccom_vcdb.npy')#np.zeros((length,100)) 
sigmat1 = np.sqrt(sigma1**2 - sigmar1**2) 
sigmat2 = np.sqrt(sigma2**2 - sigmar2**2) 
sigmat3 = np.sqrt(sigma3**2 - sigmar3**2) 
sigmat4 = np.sqrt(sigma4**2 - sigmar4**2) 

def vdisps(i, rmin=1, rmax=1e3): #sigma1 = sigma_3d_cdb_vcom, sigmar1 = sigma_r_cdb_vcom, 
    # sigma2 = sigma_3d_cdb_vcdb, sigmar2 = sigma_r_cdb_vcdb, 
    # sigma3 = sigma_3d_ccom_vcom, sigmar3 = sigma_r_ccom_vcom, 
    # sigma4 = sigma_3d_ccom_vcdb, sigmar4 = sigma_r_ccom_vcdb, 
    time = times[i]
    h1 = halo 
    print(time)
    while h1.timestep.time_gyr > time:
            h1 = h1.previous
    print (h1.timestep.time_gyr)
    snap = h1.load()
    snap.physical_units()
    print('halo collected')

    c_db = h1['shrink_center']
    c_com = pynbody.analysis.halo.center(snap, mode='com', vel=False, retcen=True)
    vc_db = h1['Vcom']
    
    core = snap.s[pynbody.filt.Sphere(10, c_db)]
    vcom_cdb = core.mean_by_mass('vel')
    core = snap.s[pynbody.filt.Sphere(10, c_com)]
    vcom_com = core.mean_by_mass('vel')
    del(core)
    gc.collect()
    write(snap.g, i, c_db, vcom_cdb, sigma = sigma1, sigmar = sigmar1, suffix = '_cdb_vcom', rmin=rmin, rmax=rmax)
    write(snap.g, i, c_db, vc_db, sigma = sigma2, sigmar = sigmar2, suffix = '_cdb_vcdb', rmin=rmin, rmax=rmax)
    
    write(snap.g, i, c_com, vcom_com, sigma = sigma3, sigmar = sigmar3, suffix = '_ccom_vcom', rmin=rmin, rmax=rmax)
    write(snap.g, i, c_com, vc_db, sigma = sigma4, sigmar = sigmar4, suffix = '_ccom_vcdb', rmin=rmin, rmax=rmax)
    del(snap); gc.collect(); gc.collect(); gc.collect()

if __name__=="__main__":
    # sigma_3d = np.load('sigma3d_halo_com.npy')
    # sigma_r = np.load('sigma_r_halo_com.npy')

    # vdisps(0, rmin=1, rmax=1e3)
    # vdisps(1, rmin=1, rmax=1e3)
    # vdisps(2, rmin=1, rmax=1e3)
    # vdisps(3, rmin=1, rmax=1e3)
    # vdisps(4, rmin=1, rmax=1e3)
    # vdisps(5, rmin=1, rmax=1e3)
    # vdisps(6, rmin=1, rmax=1e3)
    # vdisps(7, rmin=1, rmax=1e3)
    vdisps(8, rmin=1, rmax=1e3)
    vdisps(9, rmin=1, rmax=1e3)
    vdisps(10, rmin=1, rmax=1e3)
    vdisps(11, rmin=1, rmax=1e3)
    vdisps(12, rmin=1, rmax=1e3)
    vdisps(13, rmin=1, rmax=1e3)

def plot(pair, epoch):
    fig, ax = plt.subplots(ncols = 3, sharex=True, sharey=True)
    plot_profile(sigma1[pair[0]:pair[1]], color=color(pairs[3][0]), norm=False, ax=ax[0], label = r'$c_{DB},v_{c,DB}$') 
    plot_profile(sigma2[pair[0]:pair[1]], color=color(pairs[0][0]), norm=False, ax=ax[0], label = r'$c_{DB},v_{c,COM}$')  
    plot_profile(sigma3[pair[0]:pair[1]], color=color(pairs[1][0]), norm=False, ax=ax[0], label = r'$c_{COM},v_{c,COM}$') 
    plot_profile(sigma4[pair[0]:pair[1]], color=color(pairs[2][0]), norm=False, ax=ax[0], label = r'$c_{COM},v_{c,DB}$') 
    plot_profile(sigmar1[pair[0]:pair[1]], color=color(pairs[3][0]), norm=False, ax=ax[1]) 
    plot_profile(sigmar2[pair[0]:pair[1]], color=color(pairs[0][0]), norm=False, ax=ax[1]) 
    plot_profile(sigmar3[pair[0]:pair[1]], color=color(pairs[1][0]), norm=False, ax=ax[1]) 
    plot_profile(sigmar4[pair[0]:pair[1]], color=color(pairs[2][0]), norm=False, ax=ax[1])  
    plot_profile(sigmat1[pair[0]:pair[1]], color=color(pairs[3][0]), norm=False, ax=ax[2]) 
    plot_profile(sigmat2[pair[0]:pair[1]], color=color(pairs[0][0]), norm=False, ax=ax[2]) 
    plot_profile(sigmat3[pair[0]:pair[1]], color=color(pairs[1][0]), norm=False, ax=ax[2]) 
    plot_profile(sigmat4[pair[0]:pair[1]], color=color(pairs[2][0]), norm=False, ax=ax[2])  

    h, l = ax[0].get_legend_handles_labels()
    plt.legend(h,l)
    plt.xlim(1,1e3)
    plt.xscale('log')
    plt.ylim(0,500)
    ax[0].set_ylabel(r'$\sigma_{3D} (km/s)$')
    ax[1].set_ylabel(r'$\sigma_{r} (km/s)$')
    ax[2].set_ylabel(r'$\sigma_{t} (km/s)$')
    ax[1].set_title('Epoch %d' % epoch)
    plt.show(block=False)