import glob
import numpy as np
import matplotlib
if __name__=="__main__":
      matplotlib.use('Agg') #makes sure matplotlib doesn't look for a display when plotting
import matplotlib.pylab as plt
from matplotlib import cm, colors
import gc
from scipy.stats import binned_statistic

rcore = 20
allgas = True
nbins = 25

def plot(Ks, tmin, tmax, ret_cbar=False, norm=None, m=None, ind1=4,ind2=5,xmin=None, xmax=None, ymin=None, ymax=None,log=False, suffix='', qty='K', unit = r'keV $cm^2$'):
        Ks = Ks[Ks[:,1] > 1] #mass
        Ks = Ks[Ks[:,1] < 1e15]
        Ks = Ks[np.isreal(Ks[:,ind2])]
        fig, ax = plt.subplots()
        if xmin == None:
                xmin = np.nanmin(Ks[:,ind1])
                xmax = np.nanmax(Ks[:,ind1])
        if ymin == None:
                ymin = np.nanmin(Ks[:,ind2])
                ymax = np.nanmax(Ks[:,ind2])
        if log:
                binned = binned_statistic(np.log10(Ks[:,ind1]), np.log10(Ks[:,ind2]), bins = nbins, range = (xmin, xmax))
                x = 10**binned.bin_edges[:-1]
                plt.xscale('log')
                plt.yscale('log')
                plt.ylim(10**ymin, 10**ymax)
        else:
                binned = binned_statistic(Ks[:,ind1], Ks[:,ind2], bins = nbins, range = (xmin, xmax))   
                x = binned.bin_edges[:-1]
                plt.ylim(ymin, ymax)
        median = []
        min = []
        max = []
        mcells = []
        for i in xrange(nbins):
                median.append(np.nanpercentile(Ks[:,ind2][binned.binnumber == i], 50))
                min.append(np.nanpercentile(Ks[:,ind2][binned.binnumber == i], 25))
                max.append(np.nanpercentile(Ks[:,ind2][binned.binnumber == i], 75))
                mcells.append(np.nansum(Ks[:,1][binned.binnumber == i]))

        median = np.array(median)
        min = np.array(min)
        max = np.array(max)
        mcells = np.array(mcells)
        
        # plt.plot(x, x, c='g') #ax[row,col]
        # plt.fill_between(x, min, max, color=colors,alpha=0.5)
        if ret_cbar:
                norm = colors.LogNorm(vmin = np.nanmin(mcells[mcells>0]),vmax = np.nanmax(mcells))
                m = cm.ScalarMappable(norm = norm, cmap = cm.magma)
        
        if log:
                plt.bar(x, bottom=min, height=(max - min), width=10**(binned.bin_edges[1:])-x ,color=m.to_rgba(mcells))
        else:
                plt.bar(x, bottom=min, height=(max - min), width=binned.bin_edges[1:]-x ,color=m.to_rgba(mcells))

        plt.plot(x, median, c='w')
        plt.plot(x, x, c='g')
        m.set_array([])
        plt.colorbar(m)
        plt.xlabel('%s(t=%0.2f Gyr) (%s)' % (qty, tmin, unit))
        plt.ylabel('%s(t=%0.2f Gyr) (%s)' % (qty, tmax, unit))
        plt.xlim(x.min(), x.max())
        print "plotting complete"
        plt.savefig('mixing_%0.2fGyr_%dkpc%s.png' % (tmin, rcore, suffix))
        if ret_cbar:
                return norm, m


files = glob.glob('*.npy')
files.sort()

def plotall():
        norm, m = plot(np.load(files[1]), 11.3, 11.65, ind1=4,ind2=5, xmin=0,xmax=20,xmin=0,xmax=1000,ret_cbar = True)
        plot(np.load(files[2]), 11.3, 11.65, ind1=4,ind2=5, norm=norm, m=m)
        plot(np.load(files[3]), 11.97, 12.15, ind1=4,ind2=5, norm=norm, m=m)
        plot(np.load(files[3]), 13.25, 13.55, ind1=4,ind2=5, norm=norm, m=m)
        plot(np.load(files[3]), 9.2, 9.7, ind1=4,ind2=5, norm=norm, m=m)