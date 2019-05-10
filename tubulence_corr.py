import numpy as np 
import tangos as db

simname='h1.cosmo50'
sim = db.get_simulation(simname)
steps = sim.timesteps
#black holes count as halos
#get local copies of h1.cosmo50 and the database to the office desktop
haloind = 0
startstep = -1
halo = steps[startstep].halos[0]
precip = halo.reverse_property_cascade('colder_mdot_in_profile')[0]
vdisp = halo.reverse_property_cascade('v_disp_profile')[0][-12:]
mask = np.ma.masked_greater(vdisp,0).mask
valid_precip = np.ma.masked_array(precip,-mask)
valid_vdisp = np.ma.masked_array(vdisp,-mask)
r = scipy.stats.pearsonr(valid_precip.compressed(),valid_vdisp.compressed())
#when I ran the above I got:
r = 0.36