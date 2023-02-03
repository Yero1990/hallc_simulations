import numpy as np
import LT.box as B
from LT.datafile import dfile
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy.ma as ma
from scipy import optimize

import data_handling
from data_handling import *

from matplotlib import rc
rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"

thnq_bin = 35


pm_bin, pm_avg_580s1, dataXsec_580s1, dataXsec_err_580s1 = read_hallc_data(thnq=thnq_bin, pm_set=580, data_set=1)
pm_bin, pm_avg_580s2, dataXsec_580s2, dataXsec_err_580s2 = read_hallc_data(thnq=thnq_bin, pm_set=580, data_set=2)
pm_bin, pm_avg_750s1, dataXsec_750s1, dataXsec_err_750s1 = read_hallc_data(thnq=thnq_bin, pm_set=750, data_set=1)
pm_bin, pm_avg_750s2, dataXsec_750s2, dataXsec_err_750s2 = read_hallc_data(thnq=thnq_bin, pm_set=750, data_set=2)
pm_bin, pm_avg_750s3, dataXsec_750s3, dataXsec_err_750s3 = read_hallc_data(thnq=thnq_bin, pm_set=750, data_set=3)

pm_avg_cd1, f_fsiXsec_580s1_CD = read_theoretical_models("MS_CDBonn", theory="CD-Bonn", model="FSI", thnq=thnq_bin, pm_set=580, data_set=1)
pm_avg_cd2, f_fsiXsec_580s2_CD = read_theoretical_models("MS_CDBonn", theory="CD-Bonn", model="FSI", thnq=thnq_bin, pm_set=580, data_set=2)
pm_avg_cd3, f_fsiXsec_750s1_CD = read_theoretical_models("MS_CDBonn", theory="CD-Bonn", model="FSI", thnq=thnq_bin, pm_set=750, data_set=1)
pm_avg_cd4, f_fsiXsec_750s2_CD = read_theoretical_models("MS_CDBonn", theory="CD-Bonn", model="FSI", thnq=thnq_bin, pm_set=750, data_set=2)
pm_avg_cd5, f_fsiXsec_750s3_CD = read_theoretical_models("MS_CDBonn", theory="CD-Bonn", model="FSI", thnq=thnq_bin, pm_set=750, data_set=3)


#Calculate Ratio of dataXsec to model FSI Xsec
R_cd_580s1 = dataXsec_580s1 / f_fsiXsec_580s1_CD(pm_avg_580s1)
R_cd_580s1_m = ma.masked_where(R_cd_580s1 < 0, R_cd_580s1 )
R_cd_err_580s1 = dataXsec_err_580s1 / f_fsiXsec_580s1_CD(pm_avg_580s1)

R_cd_580s2 = dataXsec_580s2 / f_fsiXsec_580s2_CD(pm_avg_580s2)
R_cd_580s2_m = ma.masked_where(R_cd_580s2 < 0, R_cd_580s2 )
R_cd_err_580s2 = dataXsec_err_580s2 / f_fsiXsec_580s2_CD(pm_avg_580s2)

R_cd_750s1 = dataXsec_750s1 / f_fsiXsec_750s1_CD(pm_avg_750s1)
R_cd_750s1_m = ma.masked_where(R_cd_750s1 < 0, R_cd_750s1 )
R_cd_err_750s1 = dataXsec_err_750s1 / f_fsiXsec_750s1_CD(pm_avg_750s1)

R_cd_750s2 = dataXsec_750s2 / f_fsiXsec_750s2_CD(pm_avg_750s2)
R_cd_750s2_m = ma.masked_where(R_cd_750s2 < 0, R_cd_750s2 )
R_cd_err_750s2 = dataXsec_err_750s2 / f_fsiXsec_750s2_CD(pm_avg_750s2)

R_cd_750s3 = dataXsec_750s3 / f_fsiXsec_750s3_CD(pm_avg_750s3)
R_cd_err_750s3 = dataXsec_err_750s3 / f_fsiXsec_750s3_CD(pm_avg_750s3)


#setup subplots for plotting Xsec and ratios
plt.subplots(2,1, figsize=(10,10))
plt.subplots_adjust(top=0.93)
plt.suptitle(r'$^{2}$H$(e,e^{\prime}p)n$ Xsec and Ratios $\theta_{nq}=%d\pm5^{\circ}$'%(thnq_bin), fontsize=22)

plt.rc('xtick', labelsize=19) 
plt.rc('ytick', labelsize=19)
 
# charge to appropiate subplot
plt.subplot(2,1,1)

B.plot_exp(pm_avg_580s1, dataXsec_580s1, dataXsec_err_580s1, marker='o', color='m', logy=True, label='580 MeV/c (set1)')
B.plot_exp(pm_avg_cd1, f_fsiXsec_580s1_CD(pm_avg_cd1), marker='', linestyle='-', color='m', logy=True, label='MS CD-Bonn FSI (set1)')

B.plot_exp(pm_avg_580s2, dataXsec_580s2, dataXsec_err_580s2, marker='o', color='c', logy=True, label='580 MeV/c (set2)')
B.plot_exp(pm_avg_cd2, f_fsiXsec_580s2_CD(pm_avg_cd2), marker='', linestyle='-', color='c', logy=True, label='MS CD-Bonn FSI (set2)')

B.plot_exp(pm_avg_750s1, dataXsec_750s1, dataXsec_err_750s1, marker='o', color='b', logy=True, label='750 MeV/c (set1)')
B.plot_exp(pm_avg_cd3, f_fsiXsec_750s1_CD(pm_avg_cd3), marker='', linestyle='-', color='b', logy=True, label='MS CD-Bonn FSI (set1)')

B.plot_exp(pm_avg_750s2, dataXsec_750s2, dataXsec_err_750s2, marker='o', color='g', logy=True, label='750 MeV/c (set2)')
B.plot_exp(pm_avg_cd4, f_fsiXsec_750s2_CD(pm_avg_cd4), marker='', linestyle='-', color='g', logy=True, label='MS CD-Bonn FSI (set2)')

B.plot_exp(pm_avg_750s3, dataXsec_750s3, dataXsec_err_750s3, marker='o', color='r', logy=True, label='750 MeV/c (set3)')
B.plot_exp(pm_avg_cd5, f_fsiXsec_750s3_CD(pm_avg_cd5), marker='', linestyle='-', color='r', logy=True, label='MS CD-Bonn FSI (set3)')

B.pl.ylabel(r'Cross Sections ($\mu$b$\cdot$MeV$^{-1}$sr$^{-2}$)', fontsize=22)

B.pl.legend(loc='upper left')

#-------
plt.subplot(2,1,2)


B.plot_exp(pm_avg_580s1, R_cd_580s1, R_cd_err_580s1, marker='o', linestyle='', color='m')
B.plot_exp(pm_avg_580s2, R_cd_580s2, R_cd_err_580s2, marker='o', linestyle='', color='c')

B.plot_exp(pm_avg_750s1, R_cd_750s1, R_cd_err_750s1, marker='o', linestyle='', color='b')
B.plot_exp(pm_avg_750s2, R_cd_750s2, R_cd_err_750s2, marker='o', linestyle='', color='g')
B.plot_exp(pm_avg_750s3, R_cd_750s3, R_cd_err_750s3, marker='o', linestyle='', color='r')

B.pl.ylabel(r'Cross Section Ratio, $\sigma_{\textrm{\small data}}$/$\sigma_{\textrm{\small model}}$ ', fontsize=22)
B.pl.xlabel(r'Missing Momentum, $P_{m}$ (GeV/c)', fontsize=22)


B.pl.legend()
B.pl.ylim(0.1,15)
B.pl.show()
