import numpy as np
import matplotlib
# if __name__ == "__main__":
#     matplotlib.use("Agg")
import matplotlib.pylab as plt
import pynbody, tangos, glob, gc 
from tangos.examples import mergers

sim = tangos.get_simulation('h1.cosmo50')
steps = sim.timesteps
stepnum = -1
halonum = 0 #generalize these later
halo = steps[stepnum].halos[halonum]
mtree = mergers.get_mergers_of_major_progenitor(halo)

def stellar_maps(times, width=1700, suffix=''):
        h1, h2 = mtree[2][5]
        while h1.timestep.time_gyr > times.min():
                h1 = h1.previous
                h2 = h2.previous
        h1 = h1.next
        h2 = h2.next
        for time in times:
                while h1.timestep.time_gyr < time:
                        h1 = h1.next
                        h2 = h2.next
                stars1 = h1.load().s[pynbody.filt.HighPass('tform',0)]
                stars1.physical_units()
                stars1['pos'] -= h1['shrink_center']
                img = pynbody.plot.stars.render(stars1,
                  width=1700, mag_range=[13,28], r_scale=0.5, g_scale=0.75, b_scale=1.0)
                plt.savefig('stars1_%0.2fGyr.png' % time)

                stars2 = h2.load().s[pynbody.filt.HighPass('tform',0)]
                stars2.physical_units()
                stars2['pos'] -= h1['shrink_center']
                img2 = pynbody.plot.stars.render(stars2, clear=False,
                  width=1700, mag_range=[13,28], r_scale=0.5, g_scale=0.75, b_scale=1.0)
                plt.savefig('stars1_2_%0.2fGyr.png' % time)
        

