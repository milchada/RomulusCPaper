import glob
import numpy as np
import matplotlib
if __name__=="__main__":
      matplotlib.use('Agg') #makes sure matplotlib doesn't look for a display when plotting
import matplotlib.pylab as plt
from matplotlib import cm, colors
import gc
from scipy.stats import binned_statistic


def plot(file, rcore, ind1=4,ind2=5, ret_cbar=False, m=None, xmin=None, xmax=None, 
    ymin=None, ymax=None,log=False, qty='K', unit = r'keV $cm^2$', nbins = 25, suffix=''):
    """
    file:
    Name of .npy file to load

    ind1, ind2: 
    Column numbers with the data at the initial and final times, respectively

    ret_cbar: default = False
    If true, return normalisation of colorbar and scalarMappable for quantity

    m: default = None
    cm.ScalarMappable object to convert quantities in ind1, ind2 to colours

    xmin, xmax: float, optional
    Range of x-axis. Default is range of original values

    ymin, ymax: float, optional
    Range of y-axis. Default is range of final values

    log: optional, default = True
    Whether to log the quantities before histogramming

    qty: optional, default = 'K'
    Name of quantity being plotted, to label axes

    unit: optional, default = r'kev cm$^2$'
    Unit of quantity being plotted, to label axes

    nbins: int, optional, default = 25
    Number of bins along each axis

    suffix: str, optional
    Append to end of output plot name.
    """
    tmin, tmax = file.split('_')[0].split('-') 
    Ks = np.load(file)
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
    for i in range(nbins):
        median.append(np.nanpercentile(Ks[:,ind2][binned.binnumber == i], 50))
        min.append(np.nanpercentile(Ks[:,ind2][binned.binnumber == i], 25))
        max.append(np.nanpercentile(Ks[:,ind2][binned.binnumber == i], 75))
        mcells.append(np.nansum(Ks[:,1][binned.binnumber == i]))

    median = np.array(median)
    min = np.array(min)
    max = np.array(max)
    mcells = np.array(mcells)

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
    plt.xlabel('%s(t=%s Gyr) (%s)' % (qty, tmin, unit))
    plt.ylabel('%s(t=%s Gyr) (%s)' % (qty, tmax, unit))
    plt.xlim(x.min(), x.max())
    print( "plotting complete")
    plt.savefig('mixing_%sGyr_%dkpc_%s.png' % (tmin, rcore, qty))
    if ret_cbar:
        return m

def plotall():
    files = glob.glob('*20_kpc.npy')
    files.sort()

    m_r = plot(files[-1], rcore=20, ind1=4,ind2=5, xmin=0, xmax=20,ymin=0,ymax=50,
        ret_cbar = True, qty='r', unit = r'kpc')
    m_K = plot(files[-1], rcore=20, ind1=2,ind2=3, xmin=-2,xmax=2, ymin=-2,ymax=2, 
        log=True, ret_cbar = True)
    
    for file in files[:-1]:
        plot(file, rcore=20, ind1=4,ind2=5, xmin=0, xmax=20,ymin=0, ymax=100, 
            m=m_r, qty='r', unit = r'kpc')
        plot(file, rcore=20, ind1=2,ind2=3, xmin=-2,xmax=2, ymin=-2,ymax=2, 
            log=True, m=m_K)
    
    files = glob.glob('*30_kpc.npy')
    files.sort()
    
    m_r = plot(files[-1], rcore=30, ind1=4,ind2=5, xmin=0, xmax=20,ymin=0,ymax=50,
        ret_cbar = True, qty='r', unit = r'kpc')
    m_K = plot(files[-1], rcore=30, ind1=2,ind2=3, xmin=-2,xmax=2, ymin=-2,ymax=2, 
        log=True, ret_cbar = True)
    
    for file in files[:-1]:
        plot(file, rcore=30, ind1=4,ind2=5, xmin=0, xmax=20,ymin=0, ymax=100, 
            m=m_r, qty='r', unit = r'kpc')
        plot(file, rcore=30, ind1=2,ind2=3, xmin=-2,xmax=2, ymin=-2,ymax=2, 
            log=True, m=m_K)