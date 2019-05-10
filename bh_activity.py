import tangos as db
import numpy as np
import matplotlib.pylab as plt

db=tangos.get_simulation('h1.cosmo50')
steps = db.timesteps
steptime = [step.time_gyr for step in steps]
stepnum = -1
halonum = 0 #generalize these later
halo = steps[stepnum].halos[halonum]
bh = halo['BH_central'][0]

def bh_activity(bh=bh, plot=True):
	"""
	Periodicity found using autocorrelation as shown here:
	http://qingkaikong.blogspot.com/2017/01/signal-processing-finding-periodic.html
	"""
	bh_mdot_hist_history = bh.reverse_property_cascade("BH_mdot_histogram_ave")[0]
	bh_mdot_bins = bh_mdot_hist_history[0]
	bh_mdot_hists = bh_mdot_hist_history[1:]
	histlen=np.array([len(bh_mdot_hists[i]) for i in xrange(len(bh_mdot_hists))])
	bh_mdot_hists = bh_mdot_hists[histlen==len(bh_mdot_bins)]
	bh_mdot_history = np.array([np.sum(np.multiply(bh_mdot_bins, bh_mdot_hists[i])) for i in xrange(len(bh_mdot_hists))])
	time = np.array(bh.calculate_for_progenitors("t()")[0])[0:len(bh_mdot_history)]
	tmax = time[np.argmax(bh_mdot_history[:-1])]
	print "Peak BH accretion is at t = %d Gyr" % tmax
	if plot==True:	
		plt.plot(time, bh_mdot_history)
		plt.yscale('log')
		plt.xlabel('Time (Gyr)')
		plt.ylabel(r'$\dot{M}_{BH}$ (M$_\odot$ Gyr$^{-1})$')
		plt.savefig("bh_activity_%d.png" % halonum)
	return np.argmin(abs(np.array(steptime)-tmax)) #step number of peak activity
