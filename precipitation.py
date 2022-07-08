import pynbody, gc, os
import numpy as np
import tangos as db

simname='h1.cosmo50'
sim = db.get_simulation(simname)
steps = sim.timesteps

haloind = 0
startstep = -1
h = steps[startstep].halos[0]
zs = np.array(h.calculate_for_progenitors("z()")[0])
ts = np.array(h.calculate_for_progenitors("t()")[0])
#total mass is in Msol


ind01 = np.argmin(abs(ts - 12.50)) #disruption end
ind02 = np.argmin(abs(ts - 12)) # disruption start
ind03 = np.argmin(abs(ts - 11.1)) #1st pericenter passage done
ind04 = np.argmin(abs(ts - 10.8)) #merger start
ind05 = np.argmin(abs(ts - 9.7))#post quenching end
ind06 = np.argmin(abs(ts - 9.2)) #post quenching start
ind07 = np.argmin(abs(ts - 7.7)) #pre-quenching end
ind08 = np.argmin(abs(ts - 7.3)) #pre-quenching start

epochs = list(range(ind01, ind02+1)) + list(range(ind03, ind04+1)) + list(range(ind05, ind06+1)) + list(range(ind07, ind08+1))
snaps = [str(steps[-epoch+1].filename) for epoch in epochs]

def entropy(ptcls, mw=False):
    T_kev = ptcls['temp'].in_units('K')*pynbody.units.k.in_units('keV K**-1')
    n_cc = ptcls.g['rho'].in_units('m_p cm**-3')
    entropy = T_kev*pow(n_cc,-2./3)
    return entropy #*weight

def precipitation(pg, std=False):
    vdisp = np.sqrt(pg['vx_disp']**2 + pg['vy_disp']**2 + pg['vz_disp']**2).in_units('km s**-1')
    out = np.zeros(pg.nbins)
    drbins = pg['bin_edges'][1:] - pg['bin_edges'][:-1]
    for i in range(pg.nbins):
        gas = pg.sim[pg.binind[i]]
        max_entropy = max(np.nanmean(gas['entropy']) - np.nanstd(gas['entropy']), 1.)# -= np.nanstd(gas['entropy']) - subtracting std leaves no gas in center..
        cold_gas_filt = pynbody.filt.LowPass('entropy', max_entropy) #& pynbody.filt.FamilyFilter(pynbody.family.gas) & pynbody.filt.LowPass('vr', 0)
        cold_gas = gas[cold_gas_filt]
        cold_gas = cold_gas[cold_gas['vr'] < vdisp[i]]
        mdot = cold_gas['mass'] * cold_gas['vr'] / (pynbody.units.Unit("kpc")*drbins[i])
        out[i] = sum(mdot.in_units('Msol yr**-1'))
    print( "precipitation computed")
    return out

def precip_array(snap):
    halo_ptcls = pynbody.load(snap).halos(dosort=True).load_copy(1)
    halo_ptcls.physical_units()
    # print "ptcls collected"
    pynbody.analysis.halo.center(halo_ptcls)
    halo_ptcls = halo_ptcls.g[pynbody.filt.LowPass('amiga.grp', 2)]
    halo_ptcls.g['entropy'] = entropy(halo_ptcls.g)
    print( "halo centred")
    pg = pynbody.analysis.profile.Profile(halo_ptcls.g,xmin=1,xmax=1e3,type='log', nbins=100)
    print( "profile generated")
    return precipitation(pg, std=True)

if __name__ == "__main__":
    precip = np.zeros([len(snaps), 100])
    for snap in snaps:
        precip[snaps.index(snap)] = precip_array(snap)
        np.save('precip', precip)