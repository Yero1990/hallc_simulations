import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# function to find a word in a txt file and return the line the word is in
def find(word, file_path):
    with open(file_path, 'r') as file:
        # read all content of a file
        lines = file.readlines()

        # check if string present in a file
        for line in lines:
            if word in line:
                return line
            else:
                print('string does not exist in a file')
                
# Brief: script to plot 1D or 2D histograms from numerical data files
# ipython make_plots_d2po.py 200,300 4.5 H_Q2   --> passes pmiss = 200 and 300 to pm_user and Q2 = 4.5 to q2_user arrays


# define file path where numerical histogram files are stored
fpath = "yield_estimates/d2_polarized/smallFSI/histogram_data"

# --- user arguments ---
if len(sys.argv) == 1:
    pm_user = [0]
    q2_user = [0]
    hist_user = [0]
    print('No arguments passed')

else:
    pm_user = sys.argv[1].split(',')  # central missing momentum setting
    q2_user = sys.argv[2].split(',')  # central Q2 setting
    hist_user = sys.argv[3]           # user histogram (must match filename hist prefix)
    
    pm_user = [int(s) for s in pm_user] # convert string to ints
    q2_user = [float(s) for s in q2_user] # convert string to ints

    
print('pm_user:', pm_user)
print('q2_user:', q2_user)

# loop over central missing momentum setting (input by the user)
for i, ipm in enumerate(pm_user):

    # loop over central Q2 setting (input by the user)
    for j, jq2 in enumerate(q2_user):

        fname = '%s/%s_yield_d2pol_pm%s_Q2_%s.txt'%(fpath, hist_user, ipm, jq2)

        #-----------------------------
        if os.path.exists(fname):
            print('input file:',fname)
        else:
            print('The file does not exist, continue.')
            continue
        #-----------------------------
        
        df = pd.read_csv(fname, comment='#')


        # check if histo is 2D
        if "_vs_" in hist_user:
            xbins = len(df.xb[df.yb==df.yb[0]])
            ybins = len(df.yb[df.xb==df.xb[0]])
            zcont = np.array(df.zcont)
            counts =  np.sum(df.zcont)
            #xlabel = (((find('xlabel', fname)).split(':')[1]).strip()).replace('#', '\\')
            
            #hist2d = plt.hist2d(df.x0 ,df.y0, bins=xbins, ybins, weights=zcont, cmap = 'viridis', norm=LogNorm(), label='$Q^{2}=%s$ \n $P_{m}$=%s' % (jq2, ipm))
            hist2d = plt.hist2d(df.x0 ,df.y0, bins=(xbins, ybins), weights=zcont, cmap = 'viridis', norm=LogNorm() )
            plt.text(0.67*(df.x0[df.y0==df.y0[0]]).max(), 0.8*(df.y0[df.x0==df.x0[0]]).max(), '$Q^{2}=%s$ GeV$^{2}$ \n $P_{m}$=%s MeV \n (counts = %d)' % (jq2, ipm, counts), fontsize=10)
            
            plt.colorbar()
            
            
        # plot 1D histo, otherwise
        else:
            # make plot of 1D histogram
            counts =  np.sum(df.ycont)
            #df.x0, df.ycont, np.sqrt(df.ycont)
            plt.hist(df.x0, bins=len(df.xb), weights=df.ycont, alpha=0.2, density=False, label='$Q^{2}=%s$ GeV$^{2}$, $P_{m}$=%s MeV \n (counts = %d)' % (jq2, ipm, counts) )
            #plt.plot([1,2,3], [1,2,3], linestyle='None', mec='k', markersize=20, label=r'$Q^{2}=%s$, $P_{m}$=%s' % (jq2, ipm))
            plt.xlabel(r'x-label [units]')
            plt.ylabel(r'y-label [units]')
            plt.legend(frameon=False)



        
        
plt.show()
