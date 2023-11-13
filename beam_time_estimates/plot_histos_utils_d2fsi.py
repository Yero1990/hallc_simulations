import sys
import os
import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt

'''
Brief: Plotting histos utilities specialized for
d(e,e'p) fsi studies proposal
'''

def make_ratios_d2fsi(pm_set, thrq_set, plot_flag=''):

  
    if plot_flag=='ratio':
        # set figure subplots for ratio
        fig, ax = plt.subplots(6, 5, sharex='col', sharey='row')
        fig.text(0.5, 0.01, r'Recoil Angle $\theta_{rq}$ [deg]', ha='center')
        fig.text(0.01, 0.5, r'R = FSI / PWIA', va='center', rotation='vertical')
        fig.set_size_inches(12,10, forward=True)
        
    # loop over central missing momentum setting 
    for ipm in pm_set:

        # loop over central recoil angle setting for a given central momentum
        for ithrq in thrq_set: 

            histf                 = 'H_Pm_vs_thrq_yield_d2fsi_pm%d_thrq%d.txt'%(ipm, ithrq)  # histogram file with numerical info
            histos_file_path_pwia = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_pwia/%s'%(ipm, ithrq, histf)
            histos_file_path_fsi  = 'yield_estimates/d2_fsi/histogram_data/pm%d_thrq%d_fsi/%s'%(ipm, ithrq, histf)

            rel_err_thrs = 0.3   #  relative stat. error threshold for masking

            # read dataframe
            df_fsi  = pd.read_csv(histos_file_path_fsi,  comment='#')
            df_pwia = pd.read_csv(histos_file_path_pwia, comment='#')
            
            # get central bin values arrays
            thrq_bins = (df_fsi.x0).to_numpy()
            pm_bins   = (df_fsi.y0[thrq_bins==thrq_bins[0]]).to_numpy()
            
            # get bin content / bin content error
            fsi_N       = df_fsi.zcont
            fsi_Nerr    = np.sqrt(fsi_N) 
            fsi_rel_err = fsi_Nerr / fsi_N

            fsi_N    = ma.masked_where(fsi_N==0, fsi_N)
            fsi_Nerr = ma.masked_where(fsi_N==0, fsi_Nerr)
    
            pwia_N       = df_pwia.zcont
            pwia_Nerr    = np.sqrt(pwia_N)
            pwia_rel_err = pwia_Nerr / pwia_N
            
            pwia_N    = ma.masked_where(pwia_N==0, pwia_N)
            pwia_Nerr = ma.masked_where(pwia_N==0, pwia_Nerr)
    
            ratio = fsi_N / pwia_N
            ratio_err = ratio * np.sqrt((fsi_Nerr/fsi_N)**2 + (pwia_Nerr/pwia_N)**2)
            
            ratio_rel_err = ratio_err / ratio
            
            ratio = ma.masked_where(ratio_rel_err>rel_err_thrs,     ratio)
            ratio_err = ma.masked_where(ratio_rel_err>rel_err_thrs, ratio_err)

            # loop over the pm bins
            for idx, pm_bin in enumerate(pm_bins):

                if plot_flag=='ratio':
                
                    # ---- plot ratio fsi/pwia -----
                    ax = plt.subplot(6, 5, idx+1)
                    ax.errorbar(thrq_bins[df_fsi.y0==pm_bin], ratio[df_fsi.y0==pm_bin], ratio_err[df_fsi.y0==pm_bin], marker='o', linestyle='None', label=r'$\theta_{rq}=%.1f$'%ipm)
                    #if i==thrq_c[0]:
                    #    ax_ratio.set_title(r"%.0f $\pm$ %.0f MeV"%(pm_bin*1000, pm_binw*1000/2.), fontsize=10)


    plt.tight_layout()
    plt.legend()
    plt.show()

make_ratios_d2fsi([800], [28,49], 'ratio')
