import numpy as np
import numpy.ma as ma
import pandas as pd
from matplotlib.colors import LogNorm
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
from uncertainties import unumpy
import sys

# This script reads either (pwia, fsi) for a given pm setting, and overlays the
# histograms for different settings

#user arguments
print('argv[0]:', sys.argv[0])
print('argv[1]:', sys.argv[1])
print('argv[2]:', sys.argv[2])

# Examples:
# ipython make_plots_d2po.py 200 4.5        --> passes pmiss = 200 and Q2 = 4.5 to pm_user and q2_user arrays
# ipython make_plots_d2po.py 200,300 4.5    --> passes pmiss = 200 and 300 to pm_user and Q2 = 4.5 to q2_user arrays
# ipython make_plots_d2po.py 200 4.0,4.5    --> passes pmiss = 200 to pm_user and Q2 = 4.0 and 4.5 to q2_user arrays


pm_user = sys.argv[1].split(',')
q2_user = sys.argv[2].split(',')

pm_user = [int(s) for s in pm_user] # convert string to ints
q2_user = [float(s) for s in q2_user] # convert string to ints



# set the numebr of rows of cols based on the variable to be binned (thrq)
df = pd.read_csv('yield_d2pol_pm200_Q2_4.0_modelfsi_rad_0.1uA_168.0hr.txt', comment='#')
pm_bins = (df.y0)[df.x0==df.x0[0]]       # array of pm_bin centers   [GeV]
thrq_bins = (df.x0)[df.y0==df.y0[0]]  # array of thrq_bin centers [deg]
cols = round(np.sqrt((len(thrq_bins)))) + 1
rows = round(np.sqrt((len(thrq_bins))))


hist2d_flag = True
proj_flag = False

if proj_flag:
    # set figure subplots for overlay of pwia and fsi
    fig, ax = plt.subplots(nrows=rows, ncols=cols, sharex='col', sharey='row')
    fig.text(0.5, 0.01, r'Missing Momentum [GeV/c]', ha='center')
    fig.text(0.01, 0.5, r'SIMC Yield', va='center', rotation='vertical')
    fig.set_size_inches(14,8, forward=True)

        
# central missing momentum setting
#pm_c = [200, 300, 400, 500] # MeV/c
#Q2_c = ["4.0"]

pm_c = pm_user # MeV/c
Q2_c = q2_user

for i, ipm in enumerate(pm_c):
    
    for j, jq2 in enumerate(Q2_c):

        # read data file
        df_fsi = pd.read_csv('yield_d2pol_pm%i_Q2_%.1f_modelfsi_rad_0.1uA_168.0hr.txt'%(ipm, jq2), comment='#')

        
        pm_bins = (df_fsi.y0)[df_fsi.x0==df_fsi.x0[0]]       # array of pm_bin centers   [GeV]
        thrq_bins = (df_fsi.x0)[df_fsi.y0==df_fsi.y0[0]]  # array of thrq_bin centers [deg]

        if hist2d_flag:
            # plot 2d histogram
            plt.figure(i+j)
            zcont = np.array(df_fsi.zcont)
            hist2d = plt.hist2d(df_fsi.x0 ,df_fsi.y0, bins=(len(thrq_bins), len(pm_bins)), weights=zcont, cmap = 'viridis', norm=LogNorm(), label=r'$Q^{2}=%.1f$ \n $P_{m}$=%d'%(Q2_c[j],pm_c[i]))
            plt.colorbar()  
            plt.title(r'Missing Momentum vs. Recoil Angle', fontsize=15)
            plt.xlabel(r'Recoil Angle, $\theta_{rq}$ [deg]', fontsize=15)
            plt.ylabel(r'Missing Momentum, $P_{m}$ [GeV/c]', fontsize=15)
            plt.text(0.65*thrq_bins.max(), 0.8*pm_bins.max(), r"$Q^{2}=%.1f$ GeV$^{2}$"%(Q2_c[j])+"\n"+"$P_{m}$=%d MeV/c"%(pm_c[i])+"\n"+"Yield=%d"%(zcont.sum()), fontsize=15)
            plt.tight_layout()
        
        
        if proj_flag:
            
            # loop over each thrq bin
            for idx, thrq_bin in enumerate(thrq_bins):

                # extract the bin content for the corresponding bin
                counts = (df_fsi.zcont)[df_fsi.x0==thrq_bin]
                total_counts = np.round(np.sum(counts))  # integrated counts for the corresponding bin
                counts = ma.masked_where(counts==0, counts)
                counts_err = np.sqrt(counts)
                
                
                # add a new subplot iteratively using nrows and cols
                ax = plt.subplot(rows, cols, idx + 1)

                if(len(pm_user)==1 and len(q2_user)==1):
                    # --- 1D hist projections (of a user-selected Pm, Q2 setting)----
                    ax.hist(pm_bins, len(pm_bins), weights=counts, alpha=0.5, label=r'$P_{m}=%d, Q^{2}=%s$'%(pm_user[0], q2_user[0]))
                    ax.title.set_text(r'yield:%d, $\theta_{rq}:%.1f$' %(total_counts, thrq_bin))

            
                elif(len(pm_user)>1 and len(q2_user)==1):
                    # --- 1D hist projections overlayed for all Pm setting, Q2 settings
                    ax.hist(pm_bins, len(pm_bins), weights=counts, alpha=0.5, label=r'$P_{m}=%d, Q^{2}=%s$'%(ipm, jq2))
                    ax.title.set_text(r'$\theta_{rq}:%.1f$' %(thrq_bin))
                
                elif(len(pm_user)==1 and len(q2_user)>1):
                    # --- 1D hist projections overlayed for all Q2 settings of a given Pm
                    ax.hist(pm_bins, len(pm_bins), weights=counts, alpha=0.5, label=r'$P_{m}=%d, Q^{2}=%s$'%(ipm, jq2))
                    ax.title.set_text(r'$\theta_{rq}:%.1f$' %(thrq_bin))
                


if proj_flag:
    plt.show()        
    plt.tight_layout()
    plt.legend()
    plt.show()


if hist2d_flag:
    plt.show()
