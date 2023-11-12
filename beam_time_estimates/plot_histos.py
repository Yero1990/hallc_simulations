import sys
import os
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Brief: script that loops over numerical histo files
# in a specified directory, and  makes and saves the plots

def get_label(label, ifile):
    with open(ifile, "r") as fp:
        for line in fp:
            if label in line:
                label = (line.split(':')[1]).strip()
                return label
            
def make_plots(pm_user, thrq_user, model):

    # pm_user   : central missing momentum setting (e.g. 500)
    # thrq_user :  central recoil angle (e.g. 28)
    # model     : Laget 'pwia' or 'fsi' 

    histos_file_path = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_%s/'%(pm_user, thrq_user, model)

    for fname in os.listdir (histos_file_path):
        
        hist_file = histos_file_path + fname

        print('Opening file -----> ', hist_file)
        
        df = pd.read_csv(hist_file, comment='#')
        
        # check if histo is 2D
        if ("_vs_" in fname):
            print("VALID 2D HIST")
            xbins = len(df.xb[df.yb==df.yb[0]])
            ybins = len(df.yb[df.xb==df.xb[0]])
            zcont = np.array(df.zcont)
            counts =  np.sum(df.zcont)
            
            #if counts==0: break
            xlabel = get_label('xlabel', hist_file)
            ylabel = get_label('ylabel', hist_file) 
            title  = get_label('title', hist_file) 

            hist2d = plt.hist2d(df.x0 ,df.y0, bins=(xbins, ybins), weights=zcont, cmap = 'viridis', norm=mcolors.LogNorm(vmin=0.1, vmax=np.sqrt(counts) ))
            plt.xlabel(xlabel, fontsize=12)
            plt.ylabel(ylabel, fontsize=12)
            plt.title(title,   fontsize=14)
             
            plt.text(0.6*(df.x0[df.y0==df.y0[0]]).max(), 0.7*(df.y0[df.x0==df.x0[0]]).max(), r"$\theta_{rq}=%d$ deg"%(thrq_user)+"\n"+"$P_{m}$=%d MeV"%(pm_user)+"\n"+"(counts = %d)"%(counts), fontsize=12)
            plt.colorbar(extend='max')

            plt.show()

        elif ("_2Davg" in fname):
            
            print("VALID 2D AVG")
            xbins = len(df.xb[df.yb==df.yb[0]])
            ybins = len(df.yb[df.xb==df.xb[0]])
            zcont = np.array(df.zcont)

            xlabel = get_label('xlabel', hist_file)
            ylabel = get_label('ylabel', hist_file) 
            title  = get_label('title', hist_file) 

            if(zcont.min()>0):
                hist2d = plt.hist2d(df.x0 ,df.y0, weights=zcont, bins=(xbins, ybins), cmap = 'viridis', norm=mcolors.LogNorm())
            else:
                hist2d = plt.hist2d(df.x0 ,df.y0, weights=zcont, bins=(xbins, ybins), cmap = 'viridis')
                #plt.scatter(df.x0, df.y0, c=zcont, s = 1, cmap = 'viridis', vmin = zcont.min(), vmax = zcont.max())
            plt.text(0.6*(df.x0[df.y0==df.y0[0]]).max(), 0.7*(df.y0[df.x0==df.x0[0]]).max(), r"$\theta_{rq}=%d$ deg"%(thrq_user)+"\n"+"$P_{m}$=%d MeV"%(pm_user), fontsize=12)
            plt.colorbar(extend='max')
            plt.xlabel(xlabel, fontsize=12)
            plt.ylabel(ylabel, fontsize=12)
            plt.title(title,   fontsize=14)

            plt.show()


            
        elif ("_2Davg_" and "_vs_") not in fname:
            # make plot of 1D histogram
            print("VALID 1D HIST")
            counts =  np.sum(df.ycont)

            xlabel = get_label('xlabel', hist_file)
            ylabel = get_label('ylabel', hist_file) 
            title  = get_label('title', hist_file) 

      
            plt.hist(df.x0, bins=len(df.xb), weights=df.ycont, alpha=0.2, density=False, label=r"$\theta_{rq}=%d$ deg"%(thrq_user)+"\n"+"$P_{m}$=%d MeV"%(pm_user)+ "\n"+"(counts = %d)"%(counts))
            plt.xlabel(xlabel, fontsize=12)
            plt.ylabel(ylabel, fontsize=12)
            plt.title(title,   fontsize=14)
            plt.legend(frameon=False, fontsize=12)
            plt.show()

            
make_plots(800, 79, 'pwia')
