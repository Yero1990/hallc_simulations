# This utility module contains useful scripts for histogramming
# numerical data to file
#import LT.box as B
#from LT.datafile import dfile
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import pandas as pd

#Use Nice Fonts 
from matplotlib import rc
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"

def plot_1d_hist(fname_pattern='', pm_set=[], clr=[], scl_factor=[], title=r'', xlabel=r'', ylabel=r'', simc_flag=False):

    #create an empty list of filenames
    fnames = list()
    #f = [None] * len(pm_set)

    #create figure 
    fig, ax = plt.subplots(figsize=(10,8))
    
    #loop over filenames
    for i in enumerate(pm_set):

        idx = i[0]

        #append generic file name
        fnames.append(fname_pattern % (pm_set[idx]))

        #open file
        f = dfile(fnames[idx])

        #read data from file
        xb = f['xb']    # bin number
        x0 = f['x0']    # central x-bin value
        xlo = f['xlow'] # bin low edge value
        xup = f['xup']  # bin up edge value
        nbins = max(f['xb']) # total no. of bins
        xmin = min(xlo)
        xmax = max(xup)

        # the scaler factor is to scale the bin contents by beam time (hrs), if needed)
        bin_cnt = f['bin_cnt'] * scl_factor[idx]            # bin content
        bin_cnt_err = f['bin_cnt_err'] * scl_factor[idx]  # absolute error in bin content

        # if SIMC, calculate the statistical uncertainty as sqrt(bin_content)
        if(simc_flag):
            stats_err = np.sqrt(bin_cnt)
            # divide by sqrt of bin content, and mask element if result is invalid (i.e., inf or nan)
            rel_stats_err = ma.masked_invalid( np.divide(1., stats_err) )
            

        if(simc_flag):
            plt.hist(x0, bins=nbins, weights = bin_cnt , range=[xmin, xmax], color=clr[idx], edgecolor='k', alpha=0.2 )
            plt.errorbar(x0, bin_cnt, yerr=stats_err,  linestyle='none', marker='o', color=clr[idx], ecolor=clr[idx], alpha = 0.2, ms=3, label = r'%d MeV/c (%.1f hr @ 70 $\mu$A)' '\n' r'Integral = $%.2f\pm%.1f$' % (pm_set[idx], scl_factor[idx], np.sum(bin_cnt), np.sqrt(np.sum(bin_cnt))))            
        else:
            plt.hist(x0, bins=nbins, weights = bin_cnt , range=[xmin, xmax], color=clr[idx], edgecolor='k', alpha=0.2 )
            plt.errorbar(x0, bin_cnt, yerr=bin_cnt_err,  linestyle='none', marker='o', color=clr[idx], ecolor=clr[idx], alpha = 0.4, ms=3, label = r'%d MeV/c' % (pm_set[idx]))

    plt.legend(loc='upper right', fontsize=14)
    ax.set_ylim(ymin=1e-3)
    ax.set_yscale('log')
    plt.title(title, fontsize=22)
    plt.xticks(fontsize=19)
    plt.yticks(fontsize=19)
    plt.xlabel(xlabel, fontsize=19)
    plt.ylabel(ylabel, fontsize=19)
    
    plt.show()
    


def main():
    print('Main')

    #plot kinematics for the low/high pmiss settings separately
    fname_path = './histogram_data/jml_fsi_rad_70uA_1hr/'

    #4-Momentum Transfers
    #plot_1d_hist(fname_pattern=fname_path+'yield_Q2_pm%d.txt', pm_set=[120], clr=['m'], scl_factor=[1], title=r'4-Momentum Transfer $Q^{2}$', xlabel=r'$Q^{2}$ (GeV/c)$^{2}$', ylabel=r'Yield', simc_flag=True)
    #plot_1d_hist(fname_pattern=fname_path+'yield_Q2_pm%d.txt', pm_set=[700, 800, 900], clr=['b', 'g', 'r'], scl_factor=[45., 102., 204.], title=r'4-Momentum Transfer $Q^{2}$', xlabel=r'$Q^{2}$ (GeV/c)$^{2}$', ylabel=r'Yield', simc_flag=True)

    #Missing Energy
    #t  = r'Missing Energy'
    #xl = r'Missing Energy, $E_{m}$ (GeV)'
    #plot_1d_hist(fname_pattern=fname_path+'yield_Em_pm%d.txt', pm_set=[120], clr=['m'], scl_factor=[1], title = t, xlabel = xl, ylabel=r'Yield', simc_flag=True)
    #plot_1d_hist(fname_pattern=fname_path+'yield_Em_pm%d.txt', pm_set=[700, 800, 900], clr=['b', 'g', 'r'], scl_factor=[45., 102., 204.], title = t, xlabel = xl, ylabel=r'Yield', simc_flag=True)

    #Missing Momentum
    t  = r'Missing Momentum: Rates Estimates'
    xl = r'Missing Momentum, $P_{m}$ (GeV/c)'
    #plot_1d_hist(fname_pattern=fname_path+'yield_Pm_pm%d.txt', pm_set=[120], clr=['m'], scl_factor=[1], title = t, xlabel = xl, ylabel=r'Yield', simc_flag=True)
    plot_1d_hist(fname_pattern=fname_path+'yield_Pm_pm%d.txt', pm_set=[120,700,800,900], clr=['m', 'b', 'g', 'r'], scl_factor=[1., 45., 102., 204.], title = t, xlabel = xl, ylabel=r'Counts', simc_flag=True)
    
    
if __name__ == "__main__":
    main()

