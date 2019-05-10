from read_db import *
import temp_calc as tc 
import gc


Mbh_obs = 1e9*np.array([3.74,1.3,2.48,.169,1.47,4.65,3.87,3.72,.329,.137,9.09,.978,6.15,20.80,.614,.396,2.3]) #10^9 Msun
kTx_obs = [1.351,1.051,.696,.806,.473,1.070,1.329,0.835,.347,.512,3.963,.682,1.997,5.195,.598,.785,.997] #keV
M500_obs = [4.05,2.68,1.36,1.73,.72,2.76,3.95,1.83,.43,.082,23.9,1.31,7.73,37.40,1.06,1.66,2.46] #10^13 Msun

#so here I'll want to call all the halos from z=0 and select the ones with m500 in this range
#what radius is the temperature computed within?
ts=db.get_timestep('cosmo25/h1.cosmo50PLK.1536gst1bwK1BH.004096')

# def Mbh_Tx():
Mbh = np.array(ts.calculate_all('bh().BH_mass', dosort=True)) #live calculation - see doc
Mvir = np.array(ts.gather_property('Mvir', dosort=True))
massive_halos = np.where(Mvir[0]>8e12)[0]
# print massive_halos[1:]
temps = []
# haloids_with_gas = [4,12,35]

def main():
	Mbh = []
	M
	for haloind in massive_halos:
		halo_lazy = s.halos(dosort=True)[haloind] 
		if halo_lazy.properties.get("M_gas") > 0:
			rvir = halo_lazy.properties.get("Rvir") #kpc
			halo = halo_lazy.gas.load_copy()
			tx= tc.cluster_temp(halo, rvir, temp_cut=1e6) #this is in keV
			print haloind, tx
			temps.append(tx)

	plt.scatter(temps, Mbh[haloids_with_gas], c='k',label='Romulus C')
	plt.scatter(kTx_obs, Mbh_obs, c='r', label='Bogdan+ 2017')
	plt.xscale('log')
	plt.yscale('log')
	plt.legend()
	plt.xlabel(r"kT$_x$ (keV)")
	plt.ylabel(r"log M$_{BH} (M_\odot)$")
	plt.savefig("Mbh_Tx_yescore")

if __name__=='__main__':
	main()