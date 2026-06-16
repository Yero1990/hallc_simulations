import sys
import os
import numpy as np
import pandas as pd
import numpy.ma as ma
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d

plt.rcParams["font.family"] = "Times New Roman"

# Set global frame thickness (default is usually 0.8)
plt.rcParams['axes.linewidth'] = 1.8 


# 1. Set the global font family to serif and prioritize Times New Roman
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# 2. Force mathematical text (subscripts, Greek letters, minus signs) to use the same font
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'


'''
#-----------------
# Plot Reduced Cross Sections (for overlay sequentially) 
#----------------
plt.figure(figsize=(6, 7))

# Hall A deep (2011)
df_wb = pd.read_csv('pm_distributions_halla_thnq35.txt', comment='#')


# commissioning data (2018)
df_comm = pd.read_csv('pm_distributions_hallc_thnq35.txt', comment='#')


# full data (2023)
df_full = pd.read_csv('redXsec_fullexp.txt', comment='#')


df_cd = pd.read_csv('ms_cdbonn_thnq35.txt', comment='#')
df_av = pd.read_csv('ms_v18_thnq35.txt', comment='#')
df_paris = pd.read_csv('jml_paris_thnq35.txt', comment='#')
df_jvo = pd.read_csv('jvo_wjc2_thnq35.txt', comment='#')


# interpolate
f_redXsec_pwia_cd  = interp1d(df_cd.pm_avg, df_cd.red_theoryXsec_pwia,  kind='cubic', fill_value=np.nan, bounds_error=False)
f_redXsec_fsi_cd  = interp1d(df_cd.pm_avg, df_cd.red_theoryXsec_fsi,  kind='cubic', fill_value=np.nan, bounds_error=False)

f_redXsec_pwia_av18  = interp1d(df_av.pm_avg, df_av.red_theoryXsec_pwia,  kind='cubic', fill_value=np.nan, bounds_error=False)
f_redXsec_fsi_av18   = interp1d(df_av.pm_avg, df_av.red_theoryXsec_fsi,  kind='cubic', fill_value=np.nan, bounds_error=False)

f_redXsec_pwia_paris  = interp1d(df_paris.pm_avg, df_paris.red_theoryXsec_pwia,  kind='cubic', fill_value=np.nan, bounds_error=False)
f_redXsec_fsi_paris   = interp1d(df_paris.pm_avg, df_paris.red_theoryXsec_fsi,  kind='cubic', fill_value=np.nan, bounds_error=False)

f_redXsec_pwia_jvo  = interp1d(df_jvo.pm_avg, df_jvo.red_theoryXsec_pwia,  kind='cubic', fill_value=np.nan, bounds_error=False)
f_redXsec_fsi_jvo   = interp1d(df_jvo.pm_avg, df_jvo.red_theoryXsec_fsi,  kind='cubic', fill_value=np.nan, bounds_error=False)



# hall a (2011)
plt.errorbar(df_wb.pm_avg, df_wb.red_dataXsec, df_wb.red_dataXsec_err, marker='s', ms=6.5, mec='r', mfc='white', mew=2, ecolor='r', alpha=1, linestyle='', zorder=6)


# comm data (2018)
plt.errorbar(df_comm.pm_avg, df_comm.red_dataXsec, df_comm.red_dataXsec_err, marker='o', ms=6.5, color='k', alpha=1, mec='k',  linestyle='', zorder=7)

# full data (2023)
MeV2fm = 197.3**3    #convert MeV^-3 to fm^3
plt.errorbar(df_full.pm_avg, df_full.red_dataXsec_avg * MeV2fm, df_full.red_dataXsec_avg_err, ms=6.5, marker='s', color='lightgray', mew=1, alpha=1, mec='k', linestyle='', zorder=8)


# cd-bonn
plt.plot(df_cd.pm_avg, f_redXsec_pwia_cd(df_cd.pm_avg), marker='None' , color='magenta', linestyle='--', lw=2.5, label=r'', zorder=5)
plt.plot(df_cd.pm_avg, f_redXsec_fsi_cd(df_cd.pm_avg), marker='None' , color='magenta', linestyle='-',  lw=2.5, label=r'', zorder=5)

# av18
plt.plot(df_av.pm_avg, f_redXsec_pwia_av18(df_av.pm_avg), marker='None' , color='green', linestyle='--', lw=2.5, label=r'', zorder=4)
plt.plot(df_av.pm_avg, f_redXsec_fsi_av18(df_av.pm_avg), marker='None' , color='green', linestyle='-', lw=2.5, label=r'', zorder=4)

# paris
plt.plot(df_paris.pm_avg, f_redXsec_pwia_paris(df_paris.pm_avg), marker='None' , color='blue', linestyle='--', lw=2.5, label=r'', zorder=3)
plt.plot(df_paris.pm_avg, f_redXsec_fsi_paris(df_paris.pm_avg), marker='None' , color='blue', linestyle='-', lw=2.5, label=r'', zorder=3)

# jvo
plt.plot(df_jvo.pm_avg, f_redXsec_pwia_jvo(df_jvo.pm_avg), marker='None' , color='orange', linestyle='--', lw=2.5, label=r'', zorder=3)
plt.plot(df_jvo.pm_avg, f_redXsec_fsi_jvo(df_jvo.pm_avg), marker='None' , color='orange', linestyle='-', lw=2.5, label=r'', zorder=3)


plt.xlabel(r'$p_{\text{r}}$ (GeV/c)', fontsize=26)
plt.ylabel(r'$\sigma_{\text{red}}$ (fm$^{3}$)', fontsize=26)

plt.xlim(0,1.2)


plt.yscale('log')
plt.yticks([1e-7, 1e-5, 1e-3, 1e-1, 1e1])
plt.xticks([0, 0.25, 0.5, 0.75, 1])


plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)


plt.tight_layout()
plt.show()

'''

#-----------------
# Plot Angular Distributions (for overlay sequentially) 
#----------------

plt.figure(figsize=(7, 7))

pm_bin = 0.840


# define commissioning data path
histos_file_path_comm = '../../../exp_data/comm18/redXsec_HallC_pm%d_MeV.txt'%(pm_bin*1000)

# define 2023 d(e,e'p) data (from Gema Villegas)
histos_file_path_full = '../../../exp_data/full23/redXsec_fullexp_HallC_pm%d_MeV.txt'%(pm_bin*1000)


# read d(e,e'p) 2018 commissioning data
df_comm = pd.read_csv(histos_file_path_comm, comment='#')
# read d(e,e'p) 2023  full data
df_full = pd.read_csv(histos_file_path_full, comment='#')

thrq_set = [72, 60, 49]

# define a common x range (that includes full th_rq range for multiple interpolations)
theory_thrq = np.linspace(0, 90, num=200, endpoint=True)
i = -1
clr = ['tab:blue', 'tab:orange', 'tab:green'] 
        
for ithrq in thrq_set:
   # increment index for this loop
    i = i+1
            
    # Misak Sargsian VNA calculations
    theory_file_cd  = '../../../beam_time_estimates/yield_estimates/d2_fsi/theory_calculations/q4_sig_avkin_thnq_pm/csv/csec_calc_thrq%d_3_1_1_0_12.data' %(ithrq)
    theory_file_v18 = '../../../beam_time_estimates/yield_estimates/d2_fsi/theory_calculations/q4_sig_avkin_thnq_pm/csv/csec_calc_thrq%d_2_1_1_0_12.data' %(ithrq)

    # Sabine Jeschonnek (relativistic deuteron)
    theory_file_sj = '../../../beam_time_estimates/yield_estimates/d2_fsi/theory_calculations/SJ/d2_pm800_thrq%d_fsi_rad_output_avgkin.csv' %(ithrq)

    # read the files
    df_cd  = pd.read_csv(theory_file_cd, comment='#')
    df_v18 = pd.read_csv(theory_file_v18, comment='#')

    df_sj = pd.read_csv(theory_file_sj, comment='#')
    
    pm_bw = 0.02 #bin width
    pm_min = pm_bin - pm_bw
    pm_max = pm_bin + pm_bw

    thrq_cd  = ((df_cd['th_nq_mc'])[(df_cd['pr']>pm_min) & (df_cd['pr']<pm_max) ]).to_numpy()
    ratio_cd = ((df_cd['ratio'])[(df_cd['pr']>pm_min) & (df_cd['pr']<pm_max) ]).to_numpy()

    thrq_v18  = ((df_v18['th_nq_mc'])[(df_v18['pr']>pm_min) & (df_v18['pr']<pm_max) ]).to_numpy()
    ratio_v18 = ((df_v18['ratio'])[(df_v18['pr']>pm_min) & (df_v18['pr']<pm_max) ]).to_numpy()

    thrq_sj  = ((df_sj['thnq'])[(df_sj['pmiss']>pm_min) & (df_sj['pmiss']<pm_max) ]).to_numpy()
    xsec_fsi_sj  = ((df_sj['xsec_fsi'])[(df_sj['pmiss']>pm_min) & (df_sj['pmiss']<pm_max) ]).to_numpy()
    xsec_pwia_sj = ((df_sj['xsec_pwia'])[(df_sj['pmiss']>pm_min) & (df_sj['pmiss']<pm_max) ]).to_numpy()
    ratio_sj = xsec_fsi_sj / xsec_pwia_sj
                        
    # interpolate theory curves
    if len(ratio_cd)<=3:
        f_ratio_cd  = interp1d(thrq_cd,  ratio_cd,  kind='linear', fill_value=np.nan, bounds_error=False)
    if len(ratio_v18)<=3:
        f_ratio_v18 = interp1d(thrq_v18, ratio_v18, kind='linear', fill_value=np.nan, bounds_error=False)
    if len(ratio_sj)<=3:
        f_ratio_sj = interp1d(thrq_sj, ratio_sj, kind='linear', fill_value=np.nan, bounds_error=False)

    else:                            
        f_ratio_cd  = interp1d(thrq_cd,  ratio_cd,  kind='cubic', fill_value=np.nan, bounds_error=False)
        f_ratio_v18 = interp1d(thrq_v18, ratio_v18, kind='cubic', fill_value=np.nan, bounds_error=False)
        
        f_ratio_sj = interp1d(thrq_sj, ratio_sj, kind='cubic', fill_value=np.nan, bounds_error=False)
        
        plt.plot(theory_thrq , f_ratio_cd(theory_thrq),  marker='None', color=clr[i], linestyle='--', lw=3.5, label=r'', zorder=4)
        plt.plot(theory_thrq , f_ratio_v18(theory_thrq), marker='None', color=clr[i], linestyle='-',  lw=3.5, label=r'', zorder=3)
        plt.plot(theory_thrq , f_ratio_sj(theory_thrq),  marker='None', color=clr[i], linestyle=':',  lw=3.5, label=r'', zorder=2)

                        

        # plot Laget FSI/PWIA SIMC ratios (pac 53 simulations)
        #ax.errorbar(thrq_bins[df_fsi.y0==pm_bin], ratio[df_fsi.y0==pm_bin], ratio_err[df_fsi.y0==pm_bin], marker='o', color=clr[i], mec='k', linestyle='None', ms=5, label=r'$\theta_{nq}=%.1f$ deg'%ithrq, zorder=5)
        #ax.errorbar(thrq_bins[df_fsi.y0==pm_bin], ratio[df_fsi.y0==pm_bin], ratio_err[df_fsi.y0==pm_bin], marker='o', color=clr[i], mec='k', linestyle='None', ms=5, label='', zorder=5)
        
        # avoid double plotting data
        if(ithrq==49):
            # plot  comm_data_2018/JML_Paris_PWIA SIMC ratios (for comparison)
            #ax.errorbar(df_comm.thnq, df_comm.R_paris, df_comm.R_paris_err, marker='o', mec='k', mfc='white', mew = 1.5, ecolor='k', linestyle='None', ms=7, zorder=5)
            plt.errorbar(df_comm.thnq, df_comm.R_cd,    df_comm.R_cd_err,        marker='o', ms=8.5, color='k', linestyle='None',  zorder=5, label='')
            plt.errorbar(df_full.thnq_avg, df_full.R_cd,    df_full.R_cd_err,    marker='s', ms=8.5, color='lightgray', mew = 1.,  mec='k',    linestyle='',  zorder=6, label='')


plt.xlim(20,90)
                    
plt.axhline(1, linestyle='--', color='gray', linewidth=3.5, zorder=1)                    
plt.vlines(x = 70, ymin=1, ymax=100., color = 'r', linestyle = '--', linewidth=3.5, zorder=1) # reference line at 70 deg
plt.yscale('log')
plt.ylim(-0.5,100.)

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
plt.tight_layout()
plt.show()

plt.show()

