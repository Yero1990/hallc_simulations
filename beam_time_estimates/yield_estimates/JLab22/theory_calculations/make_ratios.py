import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy
from scipy.interpolate import make_interp_spline


fname_AV18_26='d2_Eb14_Pr1_thrq26_norad_output_avgkin_AV18.csv'
fname_AV18_48='d2_Eb14_Pr1_thrq48_norad_output_avgkin_AV18.csv'
fname_AV18_65='d2_Eb14_Pr1_thrq65_norad_output_avgkin_AV18.csv'

fname_CD_26='d2_Eb14_Pr1_thrq26_norad_output_avgkin_CD-Bonn.csv'
fname_CD_48='d2_Eb14_Pr1_thrq48_norad_output_avgkin_CD-Bonn.csv'
fname_CD_65='d2_Eb14_Pr1_thrq65_norad_output_avgkin_CD-Bonn.csv'

df_av18_26 = pd.read_csv(fname_AV18_26, comment='#')
df_av18_48 = pd.read_csv(fname_AV18_48, comment='#')
df_av18_65 = pd.read_csv(fname_AV18_65, comment='#')

df_cd_26 = pd.read_csv(fname_CD_26, comment='#')
df_cd_48 = pd.read_csv(fname_CD_48, comment='#')
df_cd_65 = pd.read_csv(fname_CD_65, comment='#')


Pc_arr = [1.04, 1.6] 
# loop over different central momentum
#for Pc in np.arange(0.8, 1.2, 0.04):

for Pc in Pc_arr:

    print ('Pc = %.3f'% Pc)
    
    # require pm = 1.0 GeV/c
    df_av18_26 = df_av18_26[df_av18_26['yb']==Pc]
    df_av18_48 = df_av18_48[df_av18_48['yb']==Pc]
    df_av18_65 = df_av18_65[df_av18_65['yb']==Pc]
    
    df_cd_26 = df_cd_26[df_cd_26['yb']==Pc]
    df_cd_48 = df_cd_48[df_cd_48['yb']==Pc]
    df_cd_65 = df_cd_65[df_cd_65['yb']==Pc]
    
    # make interpolation (AV18)
    
    # th_rq = 26 deg
    X_Y_Spline_av26 = make_interp_spline(df_av18_26['xb'], df_av18_26['ratio'])
    X_av26 = np.linspace((df_av18_26['xb']).min(), (df_av18_26['xb']).max(), 25)
    Y_av26 = X_Y_Spline_av26(X_av26)
    
    # th_rq = 48 deg
    X_Y_Spline_av48 = make_interp_spline(df_av18_48['xb'], df_av18_48['ratio'])
    X_av48 = np.linspace((df_av18_48['xb']).min(), (df_av18_48['xb']).max(), 25)
    Y_av48 = X_Y_Spline_av48(X_av48)
    
    # th_rq = 65 deg
    X_Y_Spline_av65 = make_interp_spline(df_av18_65['xb'], df_av18_65['ratio'])
    X_av65 = np.linspace((df_av18_65['xb']).min(), (df_av18_65['xb']).max(), 25)
    Y_av65 = X_Y_Spline_av65(X_av65)
    
    # make interpolation (CD-Bonn)
    
    X_Y_Spline_cd26 = make_interp_spline(df_cd_26['xb'], df_cd_26['ratio'])
    X_cd26 = np.linspace((df_cd_26['xb']).min(), (df_cd_26['xb']).max(), 500)
    Y_cd26 = X_Y_Spline_cd26(X_cd26)
    
    X_Y_Spline_cd48 = make_interp_spline(df_cd_48['xb'], df_cd_48['ratio'])
    X_cd48 = np.linspace((df_cd_48['xb']).min(), (df_cd_48['xb']).max(), 500)
    Y_cd48 = X_Y_Spline_cd48(X_cd48)
    
    X_Y_Spline_cd65 = make_interp_spline(df_cd_65['xb'], df_cd_65['ratio'])
    X_cd65 = np.linspace((df_cd_65['xb']).min(), (df_cd_65['xb']).max(), 500)
    Y_cd65 = X_Y_Spline_cd65(X_cd65)
    
    
    plt.plot(X_av26, Y_av26, marker='', color='g', mec='k', linestyle='solid', markersize=8, label='AV18 26 deg')
    plt.plot(X_av48, Y_av48, marker='', color='r', mec='k', linestyle='solid', markersize=8, label='AV18 48 deg')
    plt.plot(X_av65, Y_av65, marker='', color='b', mec='k', linestyle='solid', markersize=8, label='AV18 65 deg')
    
    plt.plot(X_cd26, Y_cd26, marker='', color='g', mec='k', linestyle='dashed', markersize=8, label='CD-Bonn 26 deg')
    plt.plot(X_cd48, Y_cd48, marker='', color='r', mec='k', linestyle='dashed', markersize=8, label='CD-Bonn 48 deg')
    plt.plot(X_cd65, Y_cd65, marker='', color='b', mec='k', linestyle='dashed', markersize=8, label='CD-Bonn 65 deg')

    
    plt.title('Angular Distributions, $P_{m}$ = %.2f GeV/c'%(Pc), fontsize=18)
    plt.xlabel(r'Recoil Angle, $\theta_{rq}$ [deg]', fontsize=16)
    plt.ylabel('R= $\sigma_{FSI}$ / $\sigma_{PWIA}$', fontsize=16)
    plt.axhline(y=1.0, c='k', linestyle='dashed', linewidth=2, zorder=0)
    plt.xlim(0, 130);
    
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.legend()
    plt.show()



