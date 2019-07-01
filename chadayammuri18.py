from read_db import *#pairs, profile_array, color, plot_profile, profilelist, norms
from astropy.cosmology import Planck15
from astropy import constants

Ez = Planck15.efunc(zs)


kT500s = (h.reverse_property_cascade('T500mw')[0] + h.reverse_property_cascade('T500ew')[0])/2.
P500 = ng_500*kT500s#mcdonald 2014 eq 4
K500 = kT500s*pow(ne_500,-2./3) #mcdonald 2014 eq 4
Tmw = h.reverse_property_cascade('Tmw_tcut_profile')[0]
Tmw *= constants.k_B.to('keV K**-1').value #for entropy to be in keV cm^2

rho_e_vol = np.array([rho*units.g/(units.cm**3) for rho in h.reverse_property_cascade('rho_e_vol_profile')[0]])
rho_g_vol = m_p*rho_e_vol*mu_e #g/cm^3

def entropy(T, rho_e):
	return T/pow(rho_e, 2/3.) 

Kmw = np.array([entropy(Tmw[elt].in_units('keV')*units.keV, rho_g_vol[elt]/(m_p*units.g)) for elt in range(len(Tmw))])
profilelist = {r'K/K$_{500}$':Kmw, r'$\rho/\rho_{crit}$':rho_g_vol, r'T/T$_{500}$':Tmw}#, 'precip':-precip}
norms = {r'K/K$_{500}$':K500, r'$\rho/\rho_{crit}$':rho_crit_a, r'T/T$_{500}$':kT500s}#, 'precip':None}
T500s = h.reverse_property_cascade('T500mw')[0]*units.keV
cs = np.sqrt((T500s/(mu*constants.m_p)).to('kpc**2 yr**-2'))

bh = h['BH_central'][0]
bh_mdot = bh.reverse_property_cascade('BH_mdot_ave')[0]
bh_ts = bh.calculate_for_progenitors('t()')[0]
bh_mass = bh.reverse_property_cascade('BH_mass')[0]

def Pbh():
	eps = 0.02
	#sound crossing time - no I really need adiabatic index for this
	power = eps*(bh_mdot)*units.Msun*constants.c**2/units.yr
	return power.value*power.unit.in_units('keV yr**-1')

def r_E(ind): #a la Tang & Churazov
	t = tbh[ind]
	mdot = bhmdot[ind]
	temp = Tmw/pynbody.units.k.in_units('keV K**-1') #back to K
	rho_e = h['rho_e_vol_profile'] #cm^-3
	rho_g = rho_e/mu
	Etherm_vol = 3./2 * pynbody.units.k.in_units('erg K**-1')*rho_g*units.cc.in_units('m^-3')*temp

	dr = 0.1*pynbody.units.kpc.in_units('m')#kpc
	r = np.arange(len(rho_e))/10. *pynbody.units.kpc.in_units('m') #kpc 
	dEtherm = Etherm_vol * 4 *np.pi*r**2 * dr 
	Etherm = np.cumsum(dEtherm) #units of ergs

def Pnorm(x, P0, c500, alpha, gamma): #Arnaud et al 2010 M500
    beta = 5.4905
    ap = 0.12 #eq 7
    apx = 0.1 - (ap+0.10)*pow(x/0.5,3)/(1+pow(x/0.5,3)) #eq 8
    px = P0/(pow(c500*x, gamma)*pow(1+pow(c500*x, alpha), (beta-gamma)/alpha)) #eq 11
    return px#*pow(M500*h70/3e14,0.12) #eq 13

tff = []
for row in range(len(cummass_profile)):
	r = np.arange(len(cummass_profile[row]))/10. * 3.085e21 #kpc to cm
	mtot = cummass_profile[row]*2e33 #Msun to g
	tff.append(np.sqrt(2*r**3/(constants.G.to('cm**3 g**-1 s**-2').value*mtot))) #output in s

tcool = h.reverse_property_cascade('tcool_tcut_mw_profile_primordial')[0] #in s
tcool_tff = np.array([tcool[elt]/tff[elt] for elt in range(len(tff))])
gas_mass = h.reverse_property_cascade('gas_mass_profile')[0]
vdisps = h.reverse_property_cascade('v_disp_profile')[0]


def profile_evolutions():
	fig, ax = plt.subplots()
	ax.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
	#fig 4a
	for pair in pairs:
		profiles = profile_array(Kmw, pair[0], pair[1], ynorm = K500*(Ez**(2./3)))
		plot_profile(profiles, color=color(pair[0]), norm=True, ax=ax)
	plt.legend(loc='best')
	plt.xlim(.003,2)
	plt.xlabel(r'R/R$_{500}$', fontsize=22)
	plt.ylabel(r'K/K$_{500} \times E(z)^{2/3}$', fontsize=22)
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim(1e-2, 10)
	fig.savefig('entropy_evolution.png')
	print('entropy profile done!')

	fig, ax = plt.subplots()
	ax.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
	for pair in pairs:
		profiles = profile_array(Tmw, pair[0], pair[1], ynorm = kT500s*(Ez**(-2./3)))
		plot_profile(profiles, color=color(pair[0]),ax=ax)
	plt.legend(loc='best')
	plt.xlim(.003,2)
	plt.xlabel(r'R/R$_{500}$', fontsize=22)
	plt.ylabel(r'T/T$_{500}\times E(z)^{-2/3}$', fontsize=22)
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim(.3, 3)
	fig.savefig('temperature_evolution.png')
	print('temperature profile done!')

	fig, ax = plt.subplots()
	ax.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
	for pair in pairs:
		profiles = profile_array(Pmw, pair[0], pair[1], ynorm = P500s*(Ez**(-2./3)))
		plot_profile(profiles, color=color(pair[0]),ax=ax)
	plt.legend(loc='best')
	plt.xlim(.003,2)
	plt.xlabel(r'R/R$_{500}$', fontsize=22)
	plt.ylabel(r'P/T$_{500}\times E(z)^{-2/3}$', fontsize=22)
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim(.3, 3)
	fig.savefig('pressure_evolution.png')
	print('pressure profile done!')

	fig, ax = plt.subplots()
	ax.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
	for pair in pairs:
		profiles = profile_array(rho_g_vol, pair[0], pair[1], ynorm = rho_crit_a)
		plot_profile(profiles, color=color(pair[0]), norm=True, ax=ax)
	plt.legend(loc='best')
	plt.xlim(.003,2)
	plt.xlabel(r'R/R$_{500}$', fontsize=22)
	plt.ylabel(r'$\rho/\rho_{crit}$', fontsize=22)
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim(5, 2e4)
	fig.savefig('density_evolution.png')
	print('density evolution done!')

	#fig 6a
	fig, ax = plt.subplots()
	ax.xaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
	precip = h.reverse_property_cascade('colder_mdot_in_profile')[0]
	precip[precip == 0.] = np.nan 
	for pair in pairs:
		plot_profile(-precip[pair[0]:pair[1]], color=color(pair[0]), norm=False, label=labels[pair[0]], ax=ax)
	plt.xlim(1,1000)
	plt.xlabel(r'R (kpc)', fontsize=22)
	plt.ylabel(r'$\dot{M(K < \bar{K}(r))}_{in}$ (M$\_odot yr^{-1}$)', fontsize=22)
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim(.1, 1e4)
	fig.savefig('precipitation_evolution.png')
	print('precipitation evolution done!')

def tcool_plots(ynorm='tff', xnorm = False):
	plt.clf()
	fig, ax = plt.subplots()
	ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%0.2f'))
	ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%0.2f'))
	for pair in pairs:
		if ynorm == 'tff':
			profiles = profile_array(tcool_tff, pair[0], pair[1], ynorm=1)
			plt.ylim(2,1e4)
			plt.ylabel(r'$t_{cool}/t_{ff}$', fontsize=22)
		elif ynorm == 'tH':
			profiles = profile_array(tcool, pair[0], pair[1], ynorm=ts*units.Gyr.in_units('s')) 
			plt.ylim(1e-3, 10)
			plt.ylabel(r'$t_{cool}/t_{H}$', fontsize=22)
		plot_profile(profiles, color = color(pair[0]), norm = True, label=labels[pair[0]],ax=ax)
	plt.xlabel(r'R/R$_{500}$', fontsize=22)		
	plt.legend(loc='best')
	plt.xlim(0.003, 2)
	plt.xscale('log')
	plt.yscale('log')
	fig.savefig('tcool_%s_evolution.png' % ynorm)
	print('tcool %s done!' % ynorm)

def vdisp_vbulk_plot(xnorm=True, ynorm=True, vdisp=False, vbulk=False, Pnt=True):
	if vdisp or Pnt:
		vdisp_hot = np.load('vdisp_hot_nosubs.npy')
		vdisps = np.empty([len(ts),100])
		vdisps[epochs] = vdisp_hot
		if vdisp:
			fig, ax = plt.subplots()
			ax.set_xlabel('R (kpc)', fontsize=22)
			ax.set_xscale('log')
			ax.set_ylabel(r'$\sqrt{\sigma_{v^2}}$ (km/s)', fontsize=22)
		if Pnt:
			fig2, ax2 = plt.subplots()
			Pthermraw = np.array([3*Tmw[elt].in_units('J')*units.J/(mu*constants.m_p.to('kg')) for elt in range(len(Tmw))])

	if vbulk:
		vbulk_hot = np.load('vbulk_hot_nosubs.npy')
		vbulks = np.empty([len(ts),100])
		vbulks[epochs] = vbulk_hot
		fig1, ax1 = plt.subplots()
		ax1.set_xlabel('R (kpc)', fontsize=22)
		ax1.set_xscale('log')
		ax1.set_ylabel(r'$\sqrt{v^2_{3D}}$ (km/s)', fontsize=22)

	for pair in pairs:
		if vdisp:
			plot_profile(vdisps[pair[0]:pair[1]], color=color(pair[0]), norm=False, ax=ax)
		if vbulk:
			plot_profile(vbulks[pair[0]:pair[1]], color=color(pair[0]), norm=False, ax=ax1)
		if Pnt:
			Pnontherm = vdisps[pair[0]:pair[1]]**2 *1e6 #m^2 s^-2
			Ptherm = profile_array(Pthermraw, pair[0], pair[1], ynorm = np.ones(len(Pthermraw)))
			Pnt_Ptot = Pnontherm/(Ptherm+Pnontherm)
			plot_profile(Pnt_Ptot, color=color(pair[0]), norm=True, ax=ax2)

	if Pnt:
		r_r500c = np.linspace(.003,2,100)
		r_r200m = r_r500c*.37 #r500c ~ 0.37r200m
		nelson_fit = 1 - .452*(1+np.exp(-(r_r200m/.841)**1.628))
		ax2.plot(r_r500c, nelson_fit, label='Nelson+ 2014',c='k')
		plt.legend()
		ax2.set_xlabel(r'R/R$_{500}$', fontsize=22)
		ax2.set_xscale('log')
		ax2.set_ylabel(r'$P_{non-therm}/P_{tot}$', fontsize=22)
		fig2.savefig('Pnt_fraction_evolution.png')

	if vdisp:
		fig.savefig('vdisp_evolution.png')
	if vbulk:
		fig1.savefig('vbulk_evolution.png')

	print('done!')

cummass = profile_array(gas_mass, 0, len(gas_mass), xnorm=False, ynorm=np.ones(len(gas_mass)),statistic='max')
	#max because gas mass is cumulative
rho_bin = np.empty([len(cummass), 100])
r = 10**np.linspace(0,3,100) * units.kpc.in_units('cm')
dr = r[1:] - r[:-1]
dr = np.insert(dr, 0, r[0])
for row in range(len(cummass)):
	mass = cummass[row][1:] - cummass[row][:-1]
	mass = np.insert(mass, 0, cummass[row][0])
	rho_bin[row] = mass*units.Msun.in_units('g')/(4*np.pi*r**2*dr)

def core_average(profile, core_radius_kpc=10):
	core_prop = np.empty(len(profile))
	r = 10**np.linspace(0,3,100)
	core_ind = np.argmin(abs(r - core_radius_kpc))
	for row in range(len(core_prop)):
		mass = cummass[row][1:] - cummass[row][:-1]
		mass = np.insert(mass, 0, cummass[row][0])
		core_prop[row] = np.nansum(profile[row][:core_ind]*mass[:core_ind])/np.sum(mass[:core_ind])
	return core_prop

def bcg_history(interactive=True, prop = 'sfh'):
	plt.clf()
	if prop == 'sfh':
		mvir = h.calculate_for_progenitors('Mvir')[0]
		sfh = h.calculate_for_progenitors('SFR_encl_25Myr')[0] #check key name
		sfh_enc = np.array([sfh[row][-1] for row in range(len(sfh))]) #so i did this for the entire halo
		mstar = h.calculate_for_progenitors('Mstar')[0]
		
		sfh_core = np.array([sfh[row][r500_ind(row)/10] for row in range(len(sfh))]) 
		mstar_profile = h.calculate_for_progenitors('star_mass_profile')[0]
		mstar_core = np.array([mstar_profile[row][r500_ind(row)/10] for row in range(len(sfh))])
		
		ssfr = sfh_core/mstar_core
		ymin=1e-12
		ymax=1e-8
	elif prop == 'ent':
		Tmw_all = h.reverse_property_cascade('Tmw_allgas_profile')[0]
		Tmw_all *= pynbody.units.k.in_units("keV K**-1")
		Kmw = np.array([entropy(Tmw_all[elt], rho_bin[elt]/m_p) for elt in range(len(Tmw_all))])
		core_prop = core_average(Kmw)
		ymin = 1
		ymax = 50
	elif prop == 'precip':
		precip = -1*h.reverse_property_cascade('colder_mdot_in_nosubs')[0]
		r = 10**np.linspace(0,3,100)
		ind1 = np.argmin(abs(r - 5))
		ind2 = np.argmin(abs(r - 30))
		core_prop = np.sum(precip[:,ind1:ind2], axis=1)
		ymin = 1e9
		ymax = 1e11

	fig, ax = plt.subplots()
	ax2 = ax.twinx()
	for pair in pairs:
		ax.fill([ts[pair[0]],ts[pair[0]],ts[pair[1]],ts[pair[1]]],[ymin,ymax,ymax,ymin],color=color(pair[0]), alpha=0.5)
	linecmap = cm.seismic
	if prop == 'sfh':
		ax.plot(ts, ssfr, c=linecmap(0.95))
		ax.set_ylabel(r'sSFR(<0.1$R_{500}$) (M_$\odot yr^{-1}/M_\odot$)', color=linecmap(0.95), fontsize=22)
		smooth = binned_statistic(bh_ts, bh_mdot, range=(0,14), bins=14)
		ax2.plot(smooth.bin_edges[:-1]+0.5, smooth.statistic, c=linecmap(0.15),ls='dotted',lw=2)
	elif prop in ['ent', 'precip']:
		ax.plot(ts[:len(core_prop)], core_prop, c=linecmap(0.95))
		mtree = mergers.get_mergers_of_major_progenitor(h)
		tmerge=[ts[np.argmin(abs(zs-mtree[0][merge]))] for merge in range(len(mtree[1]))]
		# ax2.scatter(tmerge, 1/mtree[1], marker='x', c='k')
		if prop == 'ent':
			ax.set_ylabel(r'K(<10kpc) (keV cm$^2$)', color=linecmap(0.95), fontsize=22)
		else:
			ax.set_ylabel(r'${M_{cold}}$(5-30kpc) (M$_\odot$)', color=linecmap(0.95), fontsize=22)
	
	ax2.plot(bh_ts[:len(bh_mdot)], bh_mdot, c=linecmap(0.15))
	ax.set_xlabel('t (Gyr)', fontsize=22)
	ax2.set_ylabel(r'$\dot{M_\bullet} (M_\odot yr^{-1})$', color=linecmap(0.15), fontsize=22)
	ax.set_yscale('log')
	ax2.set_yscale('log')
	ax.set_ylim(ymin, ymax)
	ax2.set_ylim(1e-2, 10)
	# ax2.set_ylim(1e43, 1e4)
	fig.savefig('bcg_%s.png' % prop)
	print('%s evolution plot done!' % prop)
	plt.show(block=False)
	if interactive:
		return fig, ax, ax2

def disruption_energy_ratio(frac_r500=0.1):
	bh_power = Pbh() #keV/yr
	ind_r500 = [int(frac_r500*r500_ind(row)) for row in range(len(Kmw))] #don't compute for earliest snapshot so sizes match up
	t_soundcross = np.array([0.1*ind_r500[row]/cs[row].value for row in range(len(cs))]) #kpc/(kpc/yr)
	bh_energy = bh_power[:len(t_soundcross)] * t_soundcross
	#or i could 
	K_01r500 = np.array([(np.nansum(Kmw[row][:ind_r500[row]]*gas_mass[row][:ind_r500[row]]/np.sum(gas_mass[row][:ind_r500[row]]))).value for row in range(len(ind_r500))])#keV cm^2
	delta_K = 150 - K_01r500
	V_01r500 = 4*np.pi*(np.array(ind_r500)/10.)**3 #kpc^3
	N_01r500 = np.array([sum(gas_mass[row][:ind_r500[row]])/(m_p*units.g.in_units('Msun')) for row in range(len(ind_r500))]) #number of particles
	rhs = 9./16 * delta_K * N_01r500**(5./3)/(V_01r500**(2./3) * units.kpc.in_units('cm')**2)
	return bh_energy/rhs #if this is > 1, E_bh is enough to raise K_core to 100kev cm^2

def plotall():
	profile_evolutions()
	vdisp_vbulk_plot()
	tcool_plots('tff')
	tcool_plots('tH')
	
