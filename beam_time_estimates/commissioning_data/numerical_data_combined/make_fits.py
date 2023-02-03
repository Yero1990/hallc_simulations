import numpy as np
from LT.datafile import dfile
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import rc
from scipy.interpolate import interp1d
from matplotlib.pyplot import *

#Use latex commands (e.g. \textit ot \textbf)
rc('text', usetex=True)
#Set default font to times new roman
font = {'family' : 'Times New Roman',
        'weight' : 'normal',
        'size'   : 12
}
plt.rc('font', **font)

#Set font
csfont = {'fontname':'Times New Roman'}

'''
#phenomenological fit 1 (exponential)
def fit_exp0(x, a, b, c):
    return a * np.exp(-b*x) + c 
p0_init_45 = (2, 1, 1e-14)


#phenomenological fit 1 (exponential + linear)
def fit_exp1(x, a, b, c, d):
    return a * np.exp(-b*x) + c*x + d 
p0_init_35 = (5.4, 8.6, 10.5, 0.6)
'''

#phenomenological fit 2 (piecewise function)
def fit_pw(x, a1,b1,c1,d, a2, b2, c2):
    return np.piecewise(x, [x<0.3, x>=0.3], [lambda x: a1*np.exp(-b1*x)+c1*x+d, lambda x: a2*np.exp(-b2*x)+c2])


def fit_pw45(x, a1,b1,c1,d, a2, b2, c2):
    return np.piecewise(x, [x<0.29, x>=0.29], [lambda x: a1*np.exp(-b1*x)+c1*x+d, lambda x: a2*np.exp(-b2*x)+c2])


p0_init_35 = (20, 30, 1e-2, 1, 1, 0.01, 1e-14)
p0_init_45 = (20, 30, 1e-2, 1, 1, 0.01, 1e-14)


def get_proj_errors(thnq):
    # Brief: This function reads the projected relative statistical errors for the upcoming full d(e,e'p)n experiment
    # These statistical errors will be projected onto the phenomenological fit to the current data

    #fname = 'd2_projected_errors_thnq%ddeg.txt' % thnq
    fname = 'd2_xsec-corr_projected_errors_thnq%ddeg.txt' % thnq

    thrs = 0.50 #require at least 50 % stats. error
    pm_bin = dfile(fname)['pm_bin']
    yields = dfile(fname)['Yield_xcorr']
    rel_stats_err = dfile(fname)['rel_stats_err_xcorr']
    return (pm_bin[(rel_stats_err<=thrs) & (yields>0)], rel_stats_err[(rel_stats_err<=thrs) & (yields>0)])

def get_Xsec(data_type='', thnq=0):

    #Brief: This function read either exp. or theory data files for plotting (re-producing our PRL plots)

    # set exp. data file names (Hall A only)
    if(data_type=='halla'):
        fname = 'pm_distributions_'+data_type+'_thnq%d.txt' % thnq
        print(fname)
        pm_bin       = dfile(fname)['pm_bin']
        pm_avg       = dfile(fname)['pm_avg']
        red_Xsec     = dfile(fname)['red_dataXsec']
        red_Xsec_err = dfile(fname)['red_dataXsec_err']
        return (pm_bin, pm_avg, red_Xsec, red_Xsec_err)

    # set exp. data file names (Hall C only)
    if(data_type=='hallc'):
        fname = 'pm_distributions_'+data_type+'_thnq%d.txt' % thnq
        print(fname)
        pm_bin       = dfile(fname)['pm_bin']
        pm_avg       = dfile(fname)['pm_avg']
        red_Xsec     = dfile(fname)['red_dataXsec']
        red_Xsec_err = dfile(fname)['red_dataXsec_err']
        rel_stats_err = dfile(fname)['rel_stats'] #relative statistical error only
        rel_syst_err = dfile(fname)['rel_syst'] #relative systematic error only        
        return (pm_bin, pm_avg, red_Xsec, red_Xsec_err, rel_stats_err, rel_syst_err)
    
    # set theory file names
    if(data_type=='jml_paris' or data_type=='jvo_wjc2' or data_type=='ms_cdbonn' or data_type=='ms_v18'):
        fname         = data_type + '_thnq%d.txt' % thnq        
        print(fname)
        pm_bin        = dfile(fname)['pm_bin']
        pm_avg        = dfile(fname)['pm_avg']
        red_Xsec_pwia = dfile(fname)['red_theoryXsec_pwia']
        red_Xsec_fsi  = dfile(fname)['red_theoryXsec_fsi']
        return (pm_bin, pm_avg, red_Xsec_pwia, red_Xsec_fsi)


def get_R(data_type='', thnq=0):
     #Brief: This function reads the RATIO of Xsec for plotting (re-producing our PRL plots)

    # set exp. data file names
    if(data_type=='hallc' or data_type=='halla'):
        fname = 'ratios/RATIO_pm_distributions_'+data_type+'_thnq%d.txt' % thnq
        print(fname)
        pm_bin       = dfile(fname)['pm_bin']
        pm_avg       = dfile(fname)['pm_avg']
        R            = dfile(fname)['R']
        R_err        = dfile(fname)['R_err']
        return (pm_bin, pm_avg, R, R_err)

    # set theory file names
    if(data_type=='jml_paris' or data_type=='jvo_wjc2' or data_type=='ms_cdbonn' or data_type=='ms_v18'):
        fname         = 'ratios/RATIO_' + data_type + '_thnq%d.txt' % thnq        
        print(fname)
        pm_bin        = dfile(fname)['pm_bin']
        pm_avg        = dfile(fname)['pm_avg']
        R_pwia        = dfile(fname)['R_pwia']
        R_fsi         = dfile(fname)['R_fsi']
        return (pm_bin, pm_avg, R_pwia, R_fsi)


# interpolate CD-Bonn PWIA for purposes of calculatinf the phenomenological and projected data ratio to CD-Bonn PWIA
f_red_pwiaXsec_CD35 = interp1d(get_Xsec('ms_cdbonn', 35)[1], get_Xsec('ms_cdbonn', 35)[2], fill_value='extrapolate', kind='cubic') 
f_red_pwiaXsec_CD45 = interp1d(get_Xsec('ms_cdbonn', 45)[1], get_Xsec('ms_cdbonn', 45)[2], fill_value='extrapolate', kind='cubic') 

    
#plot cross sections
plt.figure(figsize=(19,14))
plt.subplot(2,2,1)

#---------
# define pm-low/high selected values for plotting
#---------

# 35 deg
pm_low_sel35 = get_Xsec('hallc', 35)[1] < 0.25
pm_hi_sel35 = get_Xsec('hallc', 35)[1] > 0.480
pm_sel35 = (pm_low_sel35) | (pm_hi_sel35)

pm_low_sel35_proj = get_proj_errors(35)[0] < 0.25
pm_hi_sel35_proj = get_proj_errors(35)[0] > 0.480
pm_sel35_proj = (pm_low_sel35_proj) | (pm_hi_sel35_proj)

#select ONLY projected points to plot at FIXED value
pm_proj_flat35 = get_proj_errors(35)[0] > 0.960
#non-flat projected points (to select projected data that can be compared to meas. data) 
pm_proj_nflat35 = (pm_sel35_proj) & (get_proj_errors(35)[0] < 0.960)

# 45 deg
pm_low_sel45 = get_Xsec('hallc', 45)[1] < 0.25
pm_hi_sel45 = get_Xsec('hallc', 45)[1] > 0.520
pm_sel45 = (pm_low_sel45) | (pm_hi_sel45)

pm_low_sel45_proj = get_proj_errors(45)[0] < 0.25
pm_hi_sel45_proj = get_proj_errors(45)[0] > 0.520
pm_sel45_proj = (pm_low_sel45_proj) | (pm_hi_sel45_proj)


#non-flat projected points (to select projected data that can be compared to meas. data) 
pm_proj_nflat45 = (pm_sel45_proj) & (get_proj_errors(45)[0] < 1.04)

#select ONLY projected points to plot at FIXED value
pm_proj_flat45 = get_proj_errors(45)[0] > 1.04

#define colors
comm_data_clr = 'black'
halla_data_clr = 'forestgreen'
proj_data_clr = 'red'

label_fontsize=26
tick_fontsize=26
title_fontsize=26
mkr_size = 8
#---------------
# thnq = 35 deg
#---------------

plt.title(r'$\theta_{nq} = 35 \pm 5^{\circ}$', fontsize=title_fontsize)

plt.errorbar(get_Xsec('hallc', 35)[1][pm_sel35], get_Xsec('hallc', 35)[2][pm_sel35], get_Xsec('hallc', 35)[3][pm_sel35], linestyle='', ms=mkr_size, marker='o', color=comm_data_clr, mfc=comm_data_clr, label='Commissioning Data (Hall C)', zorder=4)
plt.errorbar(get_Xsec('halla', 35)[1], get_Xsec('halla', 35)[2], get_Xsec('halla', 35)[3], linestyle='', ms=mkr_size, marker='s', color=halla_data_clr, label='Hall A Data')
    
plt.plot(get_Xsec('jml_paris', 35)[1], get_Xsec('jml_paris', 35)[2], linestyle='--', color='b', label='JML Paris PWIA')
plt.plot(get_Xsec('jml_paris', 35)[1], get_Xsec('jml_paris', 35)[3], linestyle='-', color='b', label='JML Paris FSI')

plt.plot(get_Xsec('ms_cdbonn', 35)[1], get_Xsec('ms_cdbonn', 35)[2], linestyle='--', color='magenta', label='MS CD-Bonn PWIA')
plt.plot(get_Xsec('ms_cdbonn', 35)[1], get_Xsec('ms_cdbonn', 35)[3], linestyle='-', color='magenta', label='MS CD-Bonn FSI')

plt.plot(get_Xsec('ms_v18', 35)[1], get_Xsec('ms_v18', 35)[2], linestyle='--', color='g', label='MS AV18 PWIA')
plt.plot(get_Xsec('ms_v18', 35)[1], get_Xsec('ms_v18', 35)[3], linestyle='-', color='g', label='MS AV18 FSI')

plt.plot(get_Xsec('jvo_wjc2', 35)[1], get_Xsec('jvo_wjc2', 35)[2], linestyle='--', color='orange', label='JVO WJC2 PWBA')
plt.plot(get_Xsec('jvo_wjc2', 35)[1], get_Xsec('jvo_wjc2', 35)[3], linestyle='-', color='orange', label='JVO WJC2 FSI')


# ---- fit range for phenomenological fit -----
xfit_data = np.concatenate((get_Xsec('halla', 35)[1][0:11],     get_Xsec('hallc', 35)[1][0:7],  get_Xsec('hallc', 35)[1][11:24]), axis=0)
yfit_data = np.concatenate((get_Xsec('halla', 35)[2][0:11],     get_Xsec('hallc', 35)[2][0:7],  get_Xsec('hallc', 35)[2][11:24]), axis=0)
yfit_data_err = np.concatenate((get_Xsec('halla', 35)[3][0:11], get_Xsec('hallc', 35)[3][0:7],  get_Xsec('hallc', 35)[3][11:24]), axis=0)


# fit phenomennological piecewise function 
fit_parms_35, cov_fit_35 = curve_fit(f=fit_pw, xdata=xfit_data, ydata=yfit_data, sigma=yfit_data_err, p0=p0_init_35)

xlow = np.linspace(0., 0.26, 100)
xhi = np.linspace(0.4, 0.96, 100)

pheno_fit_35_low = fit_pw(xlow, *fit_parms_35)
pheno_fit_35_hi = fit_pw(xhi, *fit_parms_35)

# get projected Xsec for upcoming d(e,e'p)n experiment onto phenomenological fit
proj_data_35 = fit_pw(get_proj_errors(35)[0], *fit_parms_35)
proj_data_35_err = proj_data_35 * get_proj_errors(35)[1]

# plot phenomenological fit and projected data
#plt.plot(xlow, pheno_fit_35_low, linestyle='-', color=proj_data_clr, label='Phenomenological Fit ($f = Ae^{-bx} + c$)')
#plt.plot(xhi, pheno_fit_35_hi, linestyle='-', color=proj_data_clr)
#plt.errorbar(get_proj_errors(35)[0][pm_proj_nflat35]+0.01, proj_data_35[pm_proj_nflat35], proj_data_35_err[pm_proj_nflat35], linestyle='', ms=mkr_size, mfc='red', mec=proj_data_clr, marker='o', color=proj_data_clr, label=r'Projected Data', zorder=4)

# plot the isolated projected data
#plt.errorbar(get_proj_errors(35)[0][pm_proj_flat35]+0.01,  np.repeat(1e-6, len(get_proj_errors(35)[0]))[pm_proj_flat35], 1e-6 * get_proj_errors(35)[1][pm_proj_flat35], linestyle='', ms=mkr_size, mfc='red', mec=proj_data_clr, marker='o', color=proj_data_clr, zorder=4)

plt.xlim(-0.06, 1.26)
plt.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=tick_fontsize)

#label
plt.ylabel(r'Reduced Cross Sections, $\sigma_{\mathrm{red}}$ (fm$^{3}$)', fontsize=label_fontsize)

plt.yscale('log')

#----------------
# CHANGE SUBPLOT: plot relative errors from commissioning and projected
#----------------
plt.subplot(2,2,3)

#--------------------------
# plot errors with shaded region for better visuals
#--------------------------
#plot comm. data relative errors
plt.fill_between(get_Xsec('hallc', 35)[1][pm_low_sel35],  (-get_Xsec('hallc', 35)[3]/get_Xsec('hallc', 35)[2])[pm_low_sel35], (get_Xsec('hallc', 35)[3]/get_Xsec('hallc', 35)[2])[pm_low_sel35], color=comm_data_clr, alpha=0.2)
plt.fill_between(get_Xsec('hallc', 35)[1][pm_hi_sel35],  (-get_Xsec('hallc', 35)[3]/get_Xsec('hallc', 35)[2])[pm_hi_sel35], (get_Xsec('hallc', 35)[3]/get_Xsec('hallc', 35)[2])[pm_hi_sel35], color=comm_data_clr, alpha=0.2)

#plot projected errors
#plt.fill_between(get_proj_errors(35)[0][pm_low_sel35_proj]+0.01, -get_proj_errors(35)[1][pm_low_sel35_proj], get_proj_errors(35)[1][pm_low_sel35_proj], color=proj_data_clr, alpha=0.3)
#plt.fill_between(get_proj_errors(35)[0][pm_hi_sel35_proj]+0.01, -get_proj_errors(35)[1][pm_hi_sel35_proj], get_proj_errors(35)[1][pm_hi_sel35_proj], color=proj_data_clr, alpha=0.3)

plt.xlim(-0.06, 1.26)
plt.ylim(-0.6, 0.8)

#horizontal dashed lines to indicate +/-10 %,  and +/- 30 % relative error
plt.axhline(y=0.10, color='k', linestyle='--', linewidth=2, label=r'$\pm10\%$')
plt.axhline(y=-0.10, color='k', linestyle='--', linewidth=2)

plt.axhline(y=0.30, color='k', linestyle='-.', linewidth=2, label=r'$\pm30\%$')
plt.axhline(y=-0.30, color='k', linestyle='-.', linewidth=2)

plt.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=tick_fontsize)

# label
plt.ylabel(r'Relative Error', fontsize=label_fontsize)
plt.xlabel(r'Missing Momenta, $p_{\mathrm{m}}$ (GeV/c)', fontsize=label_fontsize)

#----------------
# thnq = 45 deg 
#----------------

plt.subplot(2,2,2)

plt.title(r'$\theta_{nq} = 45 \pm 5^{\circ}$', fontsize=title_fontsize)

plt.errorbar(get_Xsec('hallc', 45)[1][pm_sel45], get_Xsec('hallc', 45)[2][pm_sel45], get_Xsec('hallc', 45)[3][pm_sel45], linestyle='', ms=mkr_size, marker='o', color=comm_data_clr, mfc=comm_data_clr,label='Commissioning Data (Hall C)', zorder=4)
plt.errorbar(get_Xsec('halla', 45)[1], get_Xsec('halla', 45)[2], get_Xsec('halla', 45)[3], linestyle='', ms=mkr_size, marker='s', color=halla_data_clr, label='Hall A Data')

plt.plot(get_Xsec('jml_paris', 45)[1], get_Xsec('jml_paris', 45)[2], linestyle='--', color='b', label='JML Paris PWIA')
plt.plot(get_Xsec('jml_paris', 45)[1], get_Xsec('jml_paris', 45)[3], linestyle='-', color='b', label='JML Paris FSI')

plt.plot(get_Xsec('ms_cdbonn', 45)[1], get_Xsec('ms_cdbonn', 45)[2], linestyle='--', color='magenta', label='MS CD-Bonn PWIA')
plt.plot(get_Xsec('ms_cdbonn', 45)[1], get_Xsec('ms_cdbonn', 45)[3], linestyle='-', color='magenta', label='MS CD-Bonn FSI')

plt.plot(get_Xsec('ms_v18', 45)[1], get_Xsec('ms_v18', 45)[2], linestyle='--', color='g', label='MS AV18 PWIA')
plt.plot(get_Xsec('ms_v18', 45)[1], get_Xsec('ms_v18', 45)[3], linestyle='-', color='g', label='MS AV18 FSI')

plt.plot(get_Xsec('jvo_wjc2', 45)[1], get_Xsec('jvo_wjc2', 45)[2], linestyle='--', color='orange', label='JVO WJC2 PWBA')
plt.plot(get_Xsec('jvo_wjc2', 45)[1], get_Xsec('jvo_wjc2', 45)[3], linestyle='-', color='orange', label='JVO WJC2 FSI')

# concatenate data from both halls a and c for improved phenomenological fit
#xfit_data = np.concatenate((get_Xsec('hallc', 45)[1][0:7], get_Xsec('halla', 45)[1][7:13], get_Xsec('hallc', 45)[1][12:26]), axis=0)
#yfit_data = np.concatenate((get_Xsec('hallc', 45)[2][0:7], get_Xsec('halla', 45)[2][7:13], get_Xsec('hallc', 45)[2][12:26]), axis=0)
#yfit_data_err = np.concatenate((get_Xsec('hallc', 45)[3][0:7], get_Xsec('halla', 45)[3][7:13], get_Xsec('hallc', 45)[3][12:26]), axis=0)

# -- TEST FIT ----
xfit_data45 = np.concatenate((get_Xsec('halla', 45)[1][0:11],     get_Xsec('hallc', 45)[1][0:7],  get_Xsec('hallc', 45)[1][12:26]), axis=0)
yfit_data45 = np.concatenate((get_Xsec('halla', 45)[2][0:11],     get_Xsec('hallc', 45)[2][0:7],  get_Xsec('hallc', 45)[2][12:26]), axis=0)
yfit_data_err45 = np.concatenate((get_Xsec('halla', 45)[3][0:11], get_Xsec('hallc', 45)[3][0:7],  get_Xsec('hallc', 45)[3][12:26]), axis=0)

# fit phenomennological piecewise function 
fit_parms_45, cov_fit_45 = curve_fit(f=fit_pw45, xdata=xfit_data45, ydata=yfit_data45, sigma=yfit_data_err45, p0=p0_init_45)

#x = np.linspace(0., 1.05, 90)
xlow = np.linspace(0., 0.26, 100)
xhi = np.linspace(0.4, 1.05, 100)

pheno_fit_45_low = fit_pw45(xlow, *fit_parms_45)
pheno_fit_45_hi = fit_pw45(xhi, *fit_parms_45)


# get projected Xsec for upcoming d(e,e'p)n experiment onto phenomenological fit
proj_data_45 = fit_pw45(get_proj_errors(45)[0], *fit_parms_45)
proj_data_45_err = proj_data_45 * get_proj_errors(45)[1]

# plot phenomenological fit and projected data
#plt.plot(xlow, pheno_fit_45_low, linestyle='-', color=proj_data_clr, label='Phenomenological Fit ($f = Ae^{-bx} + c$)')
#plt.plot(xhi, pheno_fit_45_hi, linestyle='-', color=proj_data_clr)
#plt.errorbar(get_proj_errors(45)[0][pm_proj_nflat45]+0.01, proj_data_45[pm_proj_nflat45], proj_data_45_err[pm_proj_nflat45], linestyle='', ms=mkr_size, mfc='red', mec=proj_data_clr, marker='o', color=proj_data_clr, label=r'Projected Data', zorder=4)

# plot the isolated projected data
#plt.errorbar(get_proj_errors(45)[0][pm_proj_flat45]+0.01,  np.repeat(1e-6, len(get_proj_errors(45)[0]))[pm_proj_flat45], 1e-6 * get_proj_errors(45)[1][pm_proj_flat45], linestyle='', ms=mkr_size, mfc='red', mec=proj_data_clr, marker='o', color=proj_data_clr, zorder=4)

plt.legend(fontsize=13.5)
plt.xlim(-0.06, 1.26)

plt.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=tick_fontsize)
gca().tick_params(axis='y',labelbottom='off')
plt.yscale('log')


#plot relative errors from commissioning and projected
plt.subplot(2,2,4)

#--------------------------
# plot errors with shaded region for better visuals
#--------------------------
plt.fill_between(get_Xsec('hallc', 45)[1][pm_low_sel45],  (-get_Xsec('hallc', 45)[3]/get_Xsec('hallc', 45)[2])[pm_low_sel45], (get_Xsec('hallc', 45)[3]/get_Xsec('hallc', 45)[2])[pm_low_sel45], color=comm_data_clr, alpha=0.2, label='Commissioning Data (total error)')
plt.fill_between(get_Xsec('hallc', 45)[1][pm_hi_sel45],  (-get_Xsec('hallc', 45)[3]/get_Xsec('hallc', 45)[2])[pm_hi_sel45], (get_Xsec('hallc', 45)[3]/get_Xsec('hallc', 45)[2])[pm_hi_sel45], color=comm_data_clr, alpha=0.2)


#plt.fill_between(get_proj_errors(45)[0][pm_low_sel45_proj]+0.01, -get_proj_errors(45)[1][pm_low_sel45_proj], get_proj_errors(45)[1][pm_low_sel45_proj], color=proj_data_clr, alpha=0.3)
#plt.fill_between(get_proj_errors(45)[0][pm_hi_sel45_proj]+0.01, -get_proj_errors(45)[1][pm_hi_sel45_proj], get_proj_errors(45)[1][pm_hi_sel45_proj], color=proj_data_clr, alpha=0.3, label=r'Projected Data (statistical error)')


plt.xlim(-0.06, 1.26)
plt.ylim(-0.6, 0.8)

#horizontal dashed lines to indicate +/-10 %, +/- 20 % and +/- 40 % relative error
plt.axhline(y=0.10, color='k', linestyle='--', linewidth=2, label=r'$\pm10\%$')
plt.axhline(y=-0.10, color='k', linestyle='--', linewidth=2)

plt.axhline(y=0.30, color='k', linestyle='-.', linewidth=2, label=r'$\pm30\%$')
plt.axhline(y=-0.30, color='k', linestyle='-.', linewidth=2)


plt.xlabel(r'Missing Momenta, $p_{\mathrm{m}}$ (GeV/c)', fontsize=label_fontsize)

plt.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=tick_fontsize)

plt.legend(fontsize=20)

plt.subplots_adjust(left=0.08, bottom=0.08, right=0.96, top=0.96)
plt.savefig('fig1.eps')
plt.savefig('fig1.pdf')

#plt.show()

#------------------
#  PLOT RATIOS
#------------------

#plot cross sections
plt.figure(figsize=(19,13))
plt.subplot(1,2,1)

#---------------
# thnq = 35 deg
#---------------

plt.errorbar(get_R('hallc', 35)[1][pm_sel35], get_R('hallc', 35)[2][pm_sel35], get_R('hallc', 35)[3][pm_sel35], linestyle='', ms=mkr_size, marker='o', color=comm_data_clr, mfc=comm_data_clr, label='Commissioning Data (Hall C)', zorder=4)
plt.errorbar(get_R('halla', 35)[1], get_R('halla', 35)[2], get_R('halla', 35)[3], linestyle='', ms=mkr_size, marker='s', color=halla_data_clr, label='Hall A Data')

plt.plot(get_R('jml_paris', 35)[1], get_R('jml_paris', 35)[2], linestyle='--', color='b', label='JML Paris PWIA')
plt.plot(get_R('jml_paris', 35)[1], get_R('jml_paris', 35)[3], linestyle='-', color='b', label='JML Paris FSI')

plt.plot(get_R('ms_cdbonn', 35)[1], get_R('ms_cdbonn', 35)[2], linestyle='--', color='magenta', label='MS CD-Bonn PWIA')
plt.plot(get_R('ms_cdbonn', 35)[1], get_R('ms_cdbonn', 35)[3], linestyle='-', color='magenta', label='MS CD-Bonn FSI')

plt.plot(get_R('ms_v18', 35)[1], get_R('ms_v18', 35)[2], linestyle='--', color='g', label='MS AV18 PWIA')
plt.plot(get_R('ms_v18', 35)[1], get_R('ms_v18', 35)[3], linestyle='-', color='g', label='MS AV18 FSI')

plt.plot(get_R('jvo_wjc2', 35)[1], get_R('jvo_wjc2', 35)[2], linestyle='--', color='orange', label='JVO WJC2 PWBA')
plt.plot(get_R('jvo_wjc2', 35)[1], get_R('jvo_wjc2', 35)[3], linestyle='-', color='orange', label='JVO WJC2 FSI')

#Calculate ratio to CD-Bonn using phenomenological fit and projected errors
R_proj_35 =  proj_data_35 / f_red_pwiaXsec_CD35(get_proj_errors(35)[0])
R_proj_35_err =  proj_data_35_err / f_red_pwiaXsec_CD35(get_proj_errors(35)[0])
R_fit_35 = fit_pw(get_Xsec('ms_cdbonn', 35)[1], *fit_parms_35) / get_Xsec('ms_cdbonn', 35)[2]

#--- plot ratio of phenomenolofical fit to cd-bon pwia

pheno_sel = np.array(get_R('ms_cdbonn', 35)[1])<=0.96

pheno_sel_low = np.array(get_R('ms_cdbonn', 35)[1]) <=0.26
pheno_sel_hi = (np.array(get_R('ms_cdbonn', 35)[1]) >=0.40) & (np.array(get_R('ms_cdbonn', 35)[1])<=0.96)

#plt.plot(get_R('ms_cdbonn', 35)[1][pheno_sel], R_fit_35[pheno_sel], linestyle='-', color=proj_data_clr, label='Phenomenological Fit')
plt.plot(get_R('ms_cdbonn', 35)[1][pheno_sel_low], R_fit_35[pheno_sel_low], linestyle='-', color=proj_data_clr, label='Phenomenological Fit')
plt.plot(get_R('ms_cdbonn', 35)[1][pheno_sel_hi], R_fit_35[pheno_sel_hi], linestyle='-', color=proj_data_clr)


#plot ratio of projected data to cd-bonn pwia
plt.errorbar(get_proj_errors(35)[0][pm_proj_nflat35]+0.01, R_proj_35[pm_proj_nflat35], R_proj_35_err[pm_proj_nflat35], linestyle='', ms=mkr_size, marker='o', color=proj_data_clr, mfc='red', label='Projected Data', zorder=5)

#plot ratio of projected data to cd-bonn pwia (for isolated projected data ONLY)
R_flat = np.repeat(5.54, len(get_proj_errors(35)[0]))[pm_proj_flat35]
R_flat_err  = R_flat * get_proj_errors(35)[1][pm_proj_flat35]
plt.errorbar(get_proj_errors(35)[0][pm_proj_flat35]+0.01, R_flat, R_flat_err , linestyle='', ms=mkr_size, marker='o', color=proj_data_clr, mfc='red', zorder=4)

plt.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=tick_fontsize)

plt.title(r'$\theta_{nq} = 35 \pm 5^{\circ}$', fontsize=title_fontsize)
plt.xlabel(r'Missing Momenta, $p_{\mathrm{m}}$ (GeV/c)', fontsize=label_fontsize)
plt.ylabel(r'$R = \sigma_{\mathrm{red}} / \sigma^{\textrm{\large CD-Bonn PWIA}}_{\mathrm{red}}$', fontsize=label_fontsize)
plt.yscale('log')

#---------------
# thnq = 45 deg
#---------------
plt.subplot(1,2,2)

plt.errorbar(get_R('hallc', 45)[1][pm_sel45], get_R('hallc', 45)[2][pm_sel45], get_R('hallc', 45)[3][pm_sel45], linestyle='', ms=mkr_size, marker='o', color=comm_data_clr, mfc=comm_data_clr, label='Commissioning Data (Hall C)', zorder=4)
plt.errorbar(get_R('halla', 45)[1], get_R('halla', 45)[2], get_R('halla', 45)[3], linestyle='', ms=mkr_size, marker='s', color=halla_data_clr, label='Hall A Data')

plt.plot(get_R('jml_paris', 45)[1], get_R('jml_paris', 45)[2], linestyle='--', color='b', label='JML Paris PWIA')
plt.plot(get_R('jml_paris', 45)[1], get_R('jml_paris', 45)[3], linestyle='-', color='b', label='JML Paris FSI')

plt.plot(get_R('ms_cdbonn', 45)[1], get_R('ms_cdbonn', 45)[2], linestyle='--', color='magenta', label='MS CD-Bonn PWIA')
plt.plot(get_R('ms_cdbonn', 45)[1], get_R('ms_cdbonn', 45)[3], linestyle='-', color='magenta', label='MS CD-Bonn FSI')

plt.plot(get_R('ms_v18', 45)[1], get_R('ms_v18', 45)[2], linestyle='--', color='g', label='MS AV18 PWIA')
plt.plot(get_R('ms_v18', 45)[1], get_R('ms_v18', 45)[3], linestyle='-', color='g', label='MS AV18 FSI')

plt.plot(get_R('jvo_wjc2', 45)[1], get_R('jvo_wjc2', 45)[2], linestyle='--', color='orange', label='JVO WJC2 PWBA')
plt.plot(get_R('jvo_wjc2', 45)[1], get_R('jvo_wjc2', 45)[3], linestyle='-', color='orange', label='JVO WJC2 FSI')

#Calculate ratio to CD-Bonn using phenomenological fit and projected errors
R_proj_45 =  proj_data_45 / f_red_pwiaXsec_CD45(get_proj_errors(45)[0])
R_proj_45_err =  proj_data_45_err / f_red_pwiaXsec_CD45(get_proj_errors(45)[0])
R_fit_45 = fit_pw45(get_Xsec('ms_cdbonn', 45)[1], *fit_parms_45) / get_Xsec('ms_cdbonn', 45)[2]

#plot ratio of projected data to cd-bonn pwia
plt.errorbar(get_proj_errors(45)[0][pm_proj_nflat45]+0.01, R_proj_45[pm_proj_nflat45], R_proj_45_err[pm_proj_nflat45], linestyle='', ms=mkr_size, marker='o', color=proj_data_clr, mfc='red', label='Projected Data', zorder=5)


# --plot ratio of phenomenolofical fit to cd-bon pwia---
pheno_sel = np.array(get_R('ms_cdbonn', 45)[1])<=1.05

pheno_sel_low = np.array(get_R('ms_cdbonn', 45)[1]) <=0.26
pheno_sel_hi = (np.array(get_R('ms_cdbonn', 45)[1]) >=0.40) & (np.array(get_R('ms_cdbonn', 45)[1])<=1.05)

#plt.plot(get_R('ms_cdbonn', 45)[1][pheno_sel], R_fit_45[pheno_sel], linestyle='-', color=proj_data_clr, label='Phenomenological Fit')
plt.plot(get_R('ms_cdbonn', 45)[1][pheno_sel_low], R_fit_45[pheno_sel_low], linestyle='-', color=proj_data_clr, label='Phenomenological Fit')
plt.plot(get_R('ms_cdbonn', 45)[1][pheno_sel_hi], R_fit_45[pheno_sel_hi], linestyle='-', color=proj_data_clr)

#plot ratio of projected data to cd-bonn pwia (for isolated projected data ONLY)
R_flat = np.repeat(20.7, len(get_proj_errors(45)[0]))[pm_proj_flat45]
R_flat_err  = R_flat * get_proj_errors(45)[1][pm_proj_flat45]
plt.errorbar(get_proj_errors(45)[0][pm_proj_flat45]+0.01, R_flat, R_flat_err , linestyle='', ms=mkr_size, marker='o', color=proj_data_clr, mfc='red', zorder=4)

plt.xlim(-0.06, 1.26)
plt.ylim(0.4, 40)

plt.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True, labelsize=tick_fontsize)

plt.title(r'$\theta_{nq} = 45 \pm 5^{\circ}$', fontsize=title_fontsize)
plt.xlabel(r'Missing Momenta, $p_{\mathrm{m}}$ (GeV/c)', fontsize=label_fontsize)
#plt.ylabel(r'$R = \sigma_{\mathrm{red}} / \sigma^{\textrm{\large CD-Bonn PWIA}}_{\mathrm{red}}$', fontsize=label_fontsize)
plt.yscale('log')

plt.legend(fontsize=20)
plt.subplots_adjust(left=0.06, bottom=0.08, right=0.98, top=0.96)

plt.savefig('fig2.pdf')
plt.savefig('fig2.eps')
#plt.show()
