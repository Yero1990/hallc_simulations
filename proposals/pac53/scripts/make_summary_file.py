import LT.box as B
from LT.datafile import dfile
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy.ma as ma
import sys                                     
import os                                                                                                       
from sys import argv  
import matplotlib
from matplotlib import rc
from matplotlib import *
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)

# C.Y. Jul 3, 2025
# This code is used to write into a summary file the different recoil angle data /  simc
# for purposes of plotting (in future studies)

#matplotlib.use('Agg')

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

def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr


#Conversion factor:  1 fm = 1/ (197 MeV),   
#The reduced cross section is in MeV^-3 units
MeV2fm = 197.3**3    #convert MeV^-3 to fm^3

     

#Get Reduced Xsec Data File
fname = './redXsec_combined.txt'
f = B.get_file(fname)

#Get Bin Information (Same info for all files)                                                                                                  
i_b = B.get_data(f, 'i_b')    #2D bin number                                                                    
i_x = B.get_data(f, 'i_x')    #x (th_nq) bin number                                                                                    
i_y = B.get_data(f, 'i_y')    #y (pmiss) bin number                                        
thnq = B.get_data(f, 'xb')      #th_nq value at bin center                                                          
pm =  B.get_data(f, 'yb')      #pmiss value at bin center   
pm_avg = B.get_data(f, 'pm_avg')
#pm_int = np.linspace(0.01, 1.1, 100)#B.get_data(f, 'pm_int')

#Get Combined Final Red Xsec for all kinematics
red_dataXsec_avg     = B.get_data(f,'red_dataXsec_avg')
red_dataXsec_avg_err = B.get_data(f,'red_dataXsec_avg_err')
red_dataXsec_avg_syst_err = B.get_data(f,'red_dataXsec_avg_syst_err')
red_dataXsec_avg_tot_err = B.get_data(f,'red_dataXsec_avg_tot_err')
red_pwiaXsec_avg     = B.get_data(f,'red_pwiaXsec_avg')
red_fsiXsec_avg      = B.get_data(f,'red_fsiXsec_avg')

#Get total relative errors / (Pm, thnq) bin for plotting (probably also write to table)
kin_syst_tot = dfile(fname)['kin_syst_tot']    #kinematic systematics
norm_syst_tot = dfile(fname)['norm_syst_tot']  #normalization systematics
tot_syst_err = dfile(fname)['tot_syst_err']    #total systematics
tot_stats_err = dfile(fname)['tot_stats_err']  #total statistical
tot_err = dfile(fname)['tot_err']              #overall error


def write_output():
    

    #Read This Experiment (Hall C) Data, and require better than 50% statistics
    red_dataXsec_avg_masked = np.ma.array(red_dataXsec_avg, mask=(red_dataXsec_avg_err>0.5*red_dataXsec_avg))
    red_dataXsec_avg_masked = np.ma.filled(red_dataXsec_avg_masked.astype(float), np.nan)


    '''
    #--------------------------------------------
    # write numerical data binned in theta_nq
    #--------------------------------------------
    
    # define central recoil angles
    thnq_c = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105]

    for thnq_i in thnq_c:

        print(thnq_i)
        #Convert to sig_reduced from MeV^-3 to fm^3
        red_pwiaXsec_avg_i = red_pwiaXsec_avg[thnq==thnq_i]*MeV2fm    #Laget PWIA


        
        #Define the Data to JML PWIA Ratio  
        R_data = red_dataXsec_avg_masked[thnq==thnq_i]*MeV2fm / red_pwiaXsec_avg_i
        R_data_err = red_dataXsec_avg_tot_err[thnq==thnq_i]*MeV2fm / red_pwiaXsec_avg_i
        
        
        #------Get Averaged Missing Momentum-----
        pmiss_avg = pm_avg[thnq==thnq_i]


        fout_name = 'redXsec_HallC_thnq%d_deg.txt'%(thnq_i)
        fout = open(fout_name, 'w')
        comment1='#This datafile contains redXsec and ratio to JML FSI from Hall C Deuteron Experiment: E12-10-003\n' 
        comment2='#Units: pm_avg [GeV/c] :: redXsec [fm^3],   theta_nq = %d +\- 5 deg \n' % (thnq_i)
        header='pm_avg,data_redXsec,data_redXsec_tot_err,R,R_err\n'
        fout.write(comment1)
        fout.write(comment2)
        fout.write(header)

        for i in range(len(pm_avg[thnq==thnq_i])):
            fout.write('%.5f %.5E  %.5E  %.5f  %.5f \n' % ((pm_avg[thnq==thnq_i])[i], (red_dataXsec_avg_masked[thnq==thnq_i]*MeV2fm)[i], (red_dataXsec_avg_tot_err[thnq==thnq_i]*MeV2fm)[i],R_data,R_data_err))

        fout.close()
    '''

        
    #--------------------------------------------
    # write numerical data binned in pm
    #--------------------------------------------
    
    # define central pmiss
    pm_c = [520, 560, 600, 640, 680, 720, 760, 800, 840, 880, 920, 960]

    
    for pm_i in pm_c:

        # set limits for missing momentum bins
        pm_min = pm_i/1000. - 0.02  # central momntum +/- 20 Mev
        pm_max = pm_i/1000. + 0.02
                
        # define pm selection
        pm_sel = (pm_avg>pm_min) & (pm_avg<pm_max)

        #Convert to sig_reduced from MeV^-3 to fm^3
        red_pwiaXsec_avg_i = red_pwiaXsec_avg[pm_sel]*MeV2fm    #Laget PWIA

        
        #----- read AV18, CD-Bonn theory calculations (used in deep 2018 analysis) ----

        # arrays to be saved for a range of thnq and a single pm bin 
        red_pwiaXsec_CD = []
        red_pwiaXsec_V18 = []
        thnq_cd_arr = []
        thnq_v18_arr = []

        # loop over all different thnq files for a specific pm_i, to construct an angular distribution
        for ithnq in (5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105):

            print('ithnq = ', ithnq)

            fname_cd_pwia =  'theory_deep_2018/PWIA_models/theoryXsec_CD-BonnPWIA_thnq%d.00_combined.data' % (ithnq)
            fname_av18_pwia =  'theory_deep_2018/PWIA_models/theoryXsec_V18PWIA_thnq%d.00_combined.data' % (ithnq)

            cd_pwia = dfile(fname_cd_pwia)
            v18_pwia = dfile(fname_av18_pwia)

            pm_avg_cd = cd_pwia['pm_avg']
            pm_avg_v18 = v18_pwia['pm_avg']

            # define pm selection for both cd /  v18
            pm_sel_cd = (pm_avg_cd>pm_min) & (pm_avg_cd<pm_max)
            pm_sel_v18 = (pm_avg_v18>pm_min) & (pm_avg_v18<pm_max)

            # append theory points to array for each thnq of a given pm bin
            if(len( (cd_pwia['red_pwiaXsec_theory'])[pm_sel_cd])) > 0:  # UNITS: MeV^-3
                red_pwiaXsec_CD.append(((cd_pwia['red_pwiaXsec_theory'])[pm_sel_cd])[0])
                thnq_cd_arr.append(ithnq)

            if(len( (v18_pwia['red_pwiaXsec_theory'])[pm_sel_v18])) > 0:  # UNITS: MeV^-3
                red_pwiaXsec_V18.append(((v18_pwia['red_pwiaXsec_theory'])[pm_sel_v18])[0])
                thnq_v18_arr.append(ithnq) 

        # interpolate redXsec vs. thnq for the particular pm_i bin
        f_red_pwiaXsec_CD = interp1d(thnq_cd_arr, red_pwiaXsec_CD, fill_value='extrapolate', kind='cubic')                        #CD-Bonn (M. Sargsian calculation)
        f_red_pwiaXsec_V18 = interp1d(thnq_v18_arr, red_pwiaXsec_V18, fill_value='extrapolate', kind='cubic')                        #AV18 (M. Sargsian calculation)
        
        #x = np.linspace(5,105,100)
        #plt.plot(x, f_red_pwiaXsec_CD(x), linestyle='--', marker = 'o', label='pm_i = %d'%(pm_i))
        #plt.legend()
        #plt.show()
        #-----------------------------------------------------------------------------------

        #------ define data Averaged Missing Momentum-----
        pmiss_avg = pm_avg[pm_sel]
        thnq_c = thnq[pm_sel]

        print('pmiss_avg = ',pmiss_avg)
        print('thnq_c = ',thnq_c)
        #Define the Data to JML PWIA Ratio  
        R_data = red_dataXsec_avg_masked[pm_sel]*MeV2fm / red_pwiaXsec_avg_i
        R_data_err = red_dataXsec_avg_tot_err[pm_sel]*MeV2fm / red_pwiaXsec_avg_i

        # Define the Data to MS CD-Bonn PWIA Ratio
        R_data_CD = red_dataXsec_avg_masked[pm_sel]*MeV2fm / ( f_red_pwiaXsec_CD(thnq_c) * MeV2fm )
        R_data_CD_err = red_dataXsec_avg_tot_err[pm_sel]*MeV2fm / ( f_red_pwiaXsec_CD(thnq_c) * MeV2fm )

        # Define the Data to MS AV18 PWIA Ratio
        R_data_V18 = red_dataXsec_avg_masked[pm_sel]*MeV2fm / ( f_red_pwiaXsec_V18(thnq_c) * MeV2fm )
        R_data_V18_err = red_dataXsec_avg_tot_err[pm_sel]*MeV2fm / ( f_red_pwiaXsec_V18(thnq_c) * MeV2fm )

        
   
        
        #print('red_dataXsec_avg_masked[pm_sel]*MeV2fm = ', red_dataXsec_avg_masked[pm_sel]*MeV2fm)
        #print('red_pwiaXsec_avg_i = ', red_pwiaXsec_avg_i)
        #print('f_red_pwiaXsec_CD(thnq_c) * MeV2fm = ', f_red_pwiaXsec_CD(thnq_c) * MeV2fm)

        #print('pm-avg = ', pmiss_avg)
        #print('pm_c = ', pm[pm_sel])
        #print('thnq_c = ',thnq_c)
        #print('R_data = ', R_data)
        #print('R_data_CD = ', R_data_CD)
        
        fout_name = '../data/redXsec_HallC_pm%d_MeV.txt'%(pm_i)
        fout = open(fout_name, 'w')
        comment1='#This datafile contains exp. redXsec and ratio to (JML Paris, MS CD-Bonn, MS AV18) PWIA, from Hall C Deuteron Experiment: E12-10-003\n' 
        comment2='#Units: pm_avg [GeV/c] :: redXsec [fm^3], pm = %d +/- 20 MeV/c \n' % (pm_i)
        header='pm_avg,pm_c,thnq,data_redXsec,data_redXsec_tot_err,redXsec_paris,redXsec_cd,redXsec_v18,R_paris,R_paris_err,R_cd,R_cd_err,R_v18,R_v18_err\n'
        fout.write(comment1)
        fout.write(comment2)
        fout.write(header)

        for i in range(len(thnq_c)):
            print(i)
            fout.write('%.5f,%.3f,%d,%.5E,%.5E,%.5E,%.5E,%.5E,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f \n' % ((pmiss_avg[i], (pm[pm_sel])[i], thnq_c[i], (red_dataXsec_avg_masked[pm_sel]*MeV2fm)[i], (red_dataXsec_avg_tot_err[pm_sel]*MeV2fm)[i], red_pwiaXsec_avg_i[i], ( f_red_pwiaXsec_CD(thnq_c) * MeV2fm )[i],  ( f_red_pwiaXsec_V18(thnq_c) * MeV2fm )[i], R_data[i],R_data_err[i], R_data_CD[i], R_data_CD_err[i], R_data_V18[i], R_data_V18_err[i])))

        fout.close()
        

def main():
    print('Entering Main . . .')

    write_output()

if __name__=="__main__":
    main()


