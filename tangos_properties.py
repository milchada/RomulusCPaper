from __future__ import absolute_import
from __future__ import print_function
import pynbody
from tangos.properties.pynbody import PynbodyPropertyCalculation
from tangos.properties.pynbody.centring import centred_calculation
import numpy as np

class VDispProfile(PynbodyPropertyCalculation):
	requires_particle_data=True
    names = "v_disp_profile"

    def get_profile(self, particle_data, existing_properties):
        return pynbody.analysis.profile.Profile(particle_data.g, min=1, max=1e3, type='log', ndim=3)

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
    	pg = self.get_profile(particle_data, existing_properties)
    	return pg['v_disp'].in_units('km s**-1')

def gamma(pg):
    out = np.zeros(len(pg))
    for i in range(pg.nbins):
        subs = pg.sim[pg.binind[i]]
        gas = subs.g
        Prand=(gas['rho'].in_units('g cm**-3')*(gas['v_disp'].in_units('cm s^-1'))**2)/2.
        Ptherm=1.5*gas['rho'].in_units('g cm**-3')*gas['temp'].in_units('K')*pynbody.units.k.in_units('J K**-1')/(gas['mu']*pynbody.units.m_p.in_units('g'))
        Ptherm*= 1e7 #1e3 kg to g, 1e4 m^2 to cm^2 in J 
        Ptherm.units = 'g cm**-1 s**-2'
        Ptot = Prand+Ptherm
        mask = np.ma.masked_invalid(Ptot).mask
        pressure = np.ma.masked_array(Ptot,mask).compressed()
        rho = np.ma.masked_array(gas['rho'], mask).compressed()
        if pressure.size == 0:
            out[i] = np.nan
        else:
            out[i] = np.polyfit(np.log10(pressure), np.log10(rho),deg = 1)[0]
    return out

class PolyTropicIndex(PynbodyPropertyCalculation):
    requires_particle_data = True
    names="gamma_profile"

    def get_profile(self, particle_data, existing_properties):
        return pynbody.analysis.profile.Profile(particle_data.g, min=1, max=1e3, type='log', ndim=3)

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
        pg = self.get_profile(particle_data, existing_properties)
        return gamma(pg)

class ThermProfiles(PynbodyPropertyCalculation):
	requires_particle_data=True
    names = "temp_1e6.5K_profile", "rho_g_1e6.5K_profile"

    def get_profile(self, particle_data, existing_properties):
    	hot_particles = particle_data.g[particle_data.g['temp'] > pow(10,6.5)]
        return pynbody.analysis.profile.Profile(hot_particles, min=1, max=1e3, type='log', ndim=3)

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
    	pg = self.get_profile(particle_data, existing_properties)
    	return pg['temp'].in_units('K'), pg['density'].in_units('g cm**-3')
