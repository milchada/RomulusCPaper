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

vdisp = np.array([])
vbulk = np.array([])
for snap in snaps:#[(len(vdisp)/100):]:
	halo_ptcls = pynbody.load(snap).halos(dosort=True).load_copy(1)
	print("particles collected")
	halo_ptcls.physical_units()
	pynbody.analysis.halo.center(halo_ptcls, mode='ssc')
	mainhalo = halo_ptcls.g[pynbody.filt.LowPass('amiga.grp',2) * pynbody.filt.HighPass('temp', 1e6)]
	pg = pynbody.analysis.profile.Profile(mainhalo, min=1, max=1e3, type='log', ndim=3)
	print("profiles made")
	del(halo_ptcls, mainhalo)
	gc.collect()
	vd_3d = np.sqrt(pg['vx_disp'].in_units('km s**-1')**2 +
	pg['vy_disp'].in_units('km s**-1')**2 + 
	pg['vz_disp'].in_units('km s**-1')**2) 
	vdisp = np.append(vdisp, vd_3d)
	vbulk = np.append(vbulk, np.sqrt(pg['v2']))
	np.save('vdisp_hot_nosubs', vdisp)
	np.save('vbulk_hot_nosubs', vbulk)

	print "snap %d done" % snaps.index(snap)
