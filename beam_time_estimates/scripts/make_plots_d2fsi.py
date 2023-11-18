import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
from uncertainties import unumpy

# This script reads both (pwia, fsi) for a given th_rq setting, and overlays the
# FSI / PWIA for all the central thrq_c, binned in missing momentum (pm_bin)

#plot_flag = 'ratio'
#plot_flag = 'ratio_rel_err'
plot_flag = 'overlay'

#thrq_c = [49, 60, 72] #[28, 49, 55, 60, 66, 72, 79 ] #deg
thrq_c = [28, 49, 55] #deg

#thrq_c = [60 ] #deg

pm_c = [800] # MeV/c

# CREATE FIGURE FOR  INDIVIDUAL PLOTS
#plt.figure(1)

# set figure subplots for overlay of pwia and fsi
fig, ax = plt.subplots(5, 8, sharex='col', sharey='row')
fig.text(0.5, 0.01, r'Recoil Angle $\theta_{rq}$ [deg]', ha='center')
fig.text(0.01, 0.5, r'FSI (or PWIA)', va='center', rotation='vertical')
fig.set_size_inches(14,8, forward=True)

# set figure subplots for ratio
fig_ratio, ax_ratio = plt.subplots(5, 8, sharex='col', sharey='row')
fig_ratio.text(0.5, 0.01, r'Recoil Angle $\theta_{rq}$ [deg]', ha='center')
fig_ratio.text(0.01, 0.5, r'R = FSI / PWIA', va='center', rotation='vertical')
fig_ratio. set_size_inches(14,8, forward=True)

# scale counts by hours
scale = 1.
# loop over central recoil angles files for a given central momentum
for i in thrq_c:

    rel_err_thrs = 0.3 # mask >30 % relative error

    # read dataframe
    df_fsi = pd.read_csv('yield_estimates/d2_fsi/histogram_data/pm800_thrq%d_fsi/H_Pm_vs_thrq_yield_d2fsi_pm800_thrq%d.txt'%(i,i), comment='#')
    df_pwia = pd.read_csv('yield_estimates/d2_fsi/histogram_data/pm800_thrq%d_pwia/H_Pm_vs_thrq_yield_d2fsi_pm800_thrq%d.txt'%(i,i), comment='#')

    thrq_bin = df_fsi.x0
    pm_bins = df_fsi.y0[thrq_bin==thrq_bin[0]] # only use a set of pm_bins (avoid duplicates)

    pm_binw = (df_fsi.yup[0] - df_fsi.ylow[0])
    
    fsi_N = df_fsi.zcont * scale
    #fsi_Nerr = df_fsi.zcont_err * scale  # error based on weight
    fsi_Nerr = np.sqrt(fsi_N)     # error based on sqrt(N)

    fsi_rel_err = fsi_Nerr / fsi_N

    pwia_N = df_pwia.zcont * scale
    #pwia_Nerr = df_pwia.zcont_err * scale    # error based on weight
    pwia_Nerr = np.sqrt(pwia_N)              # error based on sqrt(N)
    pwia_rel_err = pwia_Nerr / pwia_N         

    pwia_N = ma.masked_where(pwia_N==0, pwia_N)
    pwia_Nerr = ma.masked_where(pwia_N==0, pwia_Nerr)

    fsi_N = ma.masked_where(fsi_N==0, fsi_N)
    fsi_Nerr = ma.masked_where(fsi_N==0, fsi_Nerr)
    
    
    ratio = fsi_N / pwia_N
    ratio_err = ratio * np.sqrt((fsi_Nerr/fsi_N)**2 + (pwia_Nerr/pwia_N)**2)

    ratio_rel_err = ratio_err / ratio

    ratio = ma.masked_where(ratio_rel_err>rel_err_thrs, ratio)
    ratio_err = ma.masked_where(ratio_rel_err>rel_err_thrs, ratio_err)


    for idx, pm_bin in enumerate(pm_bins):

        # ---- plot individual figures ------
        #if pm_bin == 0.84:
            
        #    plt.figure(1)
        #    plt.errorbar(thrq_bin[df_fsi.y0==pm_bin], pwia_N[df_fsi.y0==pm_bin], pwia_Nerr[df_fsi.y0==pm_bin], marker='o', linestyle='None', mfc='black', mec='black', ecolor='black', label='PWIA')
        #    plt.errorbar(thrq_bin[df_fsi.y0==pm_bin], fsi_N[df_fsi.y0==pm_bin], fsi_Nerr[df_fsi.y0==pm_bin], marker='^', linestyle='None', mfc='red', mec='red', ecolor='red', label='FSI')

        if plot_flag=='ratio':
            # ---- plot ratio fsi/pwia -----
            ax_ratio = plt.subplot(5, 8, idx+1)
            ax_ratio.errorbar(thrq_bin[df_fsi.y0==pm_bin], ratio[df_fsi.y0==pm_bin], ratio_err[df_fsi.y0==pm_bin], marker='o', linestyle='None', label=r'$\theta_{rq}=%.1f$'%i)
            if i==thrq_c[0]:
                ax_ratio.set_title(r"%.0f $\pm$ %.0f MeV"%(pm_bin*1000, pm_binw*1000/2.), fontsize=10)

        if plot_flag=='ratio_rel_err':
            # ---- plot ratio fsi/pwia -----
            y_const = np.zeros(len(ratio_rel_err[df_fsi.y0==pm_bin]))
            yerr = ratio_rel_err[df_fsi.y0==pm_bin]
            x = thrq_bin[df_fsi.y0==pm_bin]
            
            y_const_m = ma.masked_where(ratio_rel_err[df_fsi.y0==pm_bin]>rel_err_thrs, y_const)
            yerr_m = ma.masked_where(ratio_rel_err[df_fsi.y0==pm_bin]>rel_err_thrs, yerr)
            x_m = ma.masked_where(ratio_rel_err[df_fsi.y0==pm_bin]>rel_err_thrs, x)
            print('len_yc:',len(y_const_m))
            print('len_yerr:',len(yerr_m))
            print('len_x:',len(x_m))
            
            ax_ratio = plt.subplot(5, 8, idx+1)
            ax_ratio.errorbar(x_m, y_const_m, yerr_m, marker='o', linestyle='None', label=r'$\theta_{rq}=%.1f$'%i)

            plt.axhline(y = 0.15, color = 'r', linestyle = '--')
            plt.axhline(y = -0.15, color = 'r', linestyle = '--')
            
            if i==thrq_c[0]:
                ax_ratio.set_title(r"%.0f $\pm$ %.0f MeV"%(pm_bin*1000, pm_binw*1000/2.), fontsize=10)

                
        if plot_flag=='overlay':
            print('pm_c:', pm_c[0], 'ithrq:', i, 'pm_bin:', pm_bin, 'idx:',idx)
            ax = plt.subplot(5, 8, idx+1)
            ax.errorbar(thrq_bin[df_fsi.y0==pm_bin], pwia_N[df_fsi.y0==pm_bin], pwia_Nerr[df_fsi.y0==pm_bin], marker='o', linestyle='None', mfc='None', mec='black', ecolor='black', label=r'$\theta_{rq}=%.1f$'%i)
            ax.errorbar(thrq_bin[df_fsi.y0==pm_bin], fsi_N[df_fsi.y0==pm_bin], fsi_Nerr[df_fsi.y0==pm_bin], marker='^', linestyle='None', mfc='None', mec='red', ecolor='red')
            if i==thrq_c[0]:
                ax.set_title(r"%.2f $\pm$ %.2f MeV"%(pm_bin, pm_binw/2.), fontsize=10)




            
plt.tight_layout()
plt.legend()
plt.show()

  
        
