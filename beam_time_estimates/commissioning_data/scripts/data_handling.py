'''
Author: C. Yero
Date: May 29, 2021
e-mail: cyero@jlab.org

Brief: This python module contains utility functions
specialized for handling and organizing the raw data 
and theory files containing the cross sections from the
deuteron commissioning experiment (E12-10-003).
'''

from LT.datafile import dfile
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy.ma as ma

def unpack_data(thnq=0):
    '''
    This function combines matching thnq bins from different data files 
    (i.e, pm80, pm580_set1, . . .) into a single thnq data file containing
    different overlapping pm bins. The same is also done for the theoretical 
    JML Paris cross sections, which are also present in the data files
    '''

    print('Work in progress . . . ')



    
def read_hallc_data(thnq=0, pm_set=0, data_set=0):

    #Function to read E12-10-003 experimental Xsec (or red. Xsec) per data set

    fname = "./data/pm%d_laget_bc_corr_set%d.txt" % (pm_set, data_set)

    kin = dfile(fname)
    thnq_data = np.array(kin['xb'])
    pm_bin = np.array( (kin['yb'])[thnq_data==thnq] )  #[GeV]
    pm_avg = np.array( (kin['pm_avg'])[thnq_data==thnq] ) #[GeV]
    
    #Get data average Xsec (bin-centered and radiative corrected) 
    dataXsec_avg      = np.array( (kin['fsiRC_dataXsec_fsibc_corr'])[thnq_data==thnq] )  # ub MeV^-1 sr^-2
    dataXsec_avg_err  = np.array( (kin['fsiRC_dataXsec_fsibc_corr_err'])[thnq_data==thnq] )  # absolute stats. error

    #Define array-masking condition (if error is > 50 %)
    cnd = (dataXsec_avg_err > 0.5*dataXsec_avg) 
    
    #Require better than 50% statistics
    dataXsec_avg_m = ma.masked_where(cnd, dataXsec_avg)
    dataXsec_avg_m = np.ma.filled(dataXsec_avg_m.astype(float), np.nan)
    dataXsec_avg_err_m = ma.masked_where(cnd, dataXsec_avg_err )
    dataXsec_avg_err_m = np.ma.filled(dataXsec_avg_err_m.astype(float), np.nan)
                
    return pm_bin, pm_avg, dataXsec_avg_m, dataXsec_avg_err_m   #Units: GeV/c, Gev/c, ub.MeV^-1.sr^-2, ub.MeV^-1.sr^-2 

#___________________________________________________________________
def read_theoretical_models(theory_dir="", theory="", model="", thnq=0, pm_set=0, data_set=0):

    #This code read the averaged theoretical Xsec and returns arrays in pm and interpolated Xsec
    #The interpolation values can be obtained as follows: pm_avg, f_theory = read_theoretical_models(theory="", model="", thnq=0):
    #Then, the interpolated values are:  f_theory( pm ), where pm is any missing momentum value evaluated at the function, f_theory
    #theory_dir: 
    #theory: V18, CD-Bonn  |  model: PWIA, FSI

    thnq_f = "%.2f" %(thnq)
    fname = './theory/%s/theoryXsec_%s%s_thnq%s.data' % (theory_dir, theory, model, thnq_f) #this file only reads CD-Bonn/AV18. JML will be read separately.
    try:
        kin = dfile(fname)
        setting = np.array(kin['setting'])
    except:
        print(fname,' \n does not exist.')

    
    if(theory=="V18" and model=="PWIA"):
        pwiaXsec_V18 = np.array(kin['pwiaXsec_theory'][setting=='%d_set%d'%(pm_set, data_set)]) # ub * MeV^-1 *sr^-2
        pm_avg = np.array(kin['pm_avg'][setting=='%d_set%d'%(pm_set, data_set)]) # GeV/c
        #interpolate
        f_pwiaXsec_V18 = interp1d(pm_avg, pwiaXsec_V18,fill_value='extrapolate', kind='cubic')    #AV18 (M. Sargsian calculation)
        return pm_avg, f_pwiaXsec_V18

    if(theory=="V18" and model=="FSI"):
        fsiXsec_V18 = np.array(kin['fsiXsec_theory'][setting=='%d_set%d'%(pm_set, data_set)])   # ub * MeV^-1 *sr^-2
        pm_avg = np.array(kin['pm_avg'][setting=='%d_set%d'%(pm_set, data_set)]) # GeV/c
        #interpolate
        f_fsiXsec_V18  = interp1d(pm_avg, fsiXsec_V18,fill_value='extrapolate', kind='cubic')
        return pm_avg, f_fsiXsec_V18

    
    if(theory=="CD-Bonn" and model=="PWIA"):                                     
        pwiaXsec_CD_Bonn = np.array(kin['pwiaXsec_theory'][setting=='%d_set%d'%(pm_set, data_set)])
        pm_avg = np.array(kin['pm_avg'][setting=='%d_set%d'%(pm_set, data_set)])
        #interpolate
        f_pwiaXsec_CD = interp1d(pm_avg, pwiaXsec_CD_Bonn,fill_value='extrapolate', kind='cubic') 
        return pm_avg, f_pwiaXsec_CD

    if(theory=="CD-Bonn" and model=="FSI"):                                     
        fsiXsec_CD_Bonn = np.array(kin['fsiXsec_theory'][setting=='%d_set%d'%(pm_set, data_set)])
        pm_avg = np.array(kin['pm_avg'][setting=='%d_set%d'%(pm_set, data_set)])
        #interpolate
        f_fsiXsec_CD = interp1d(pm_avg,  fsiXsec_CD_Bonn,fill_value='extrapolate', kind='cubic')        
        return pm_avg, f_fsiXsec_CD

        
    #---Reading J.M. Laget theory from the data file---
    if(theory=="JML"):

        fname = "./data/pm%d_laget_bc_corr_set%d.txt" % (pm_set, data_set)
        print(fname)
        
        kin = dfile(fname)
        thnq_data = np.array(kin['xb'])
        pm_bin = np.array( (kin['yb'])[thnq_data==thnq] )
        pm_avg = np.array( (kin['pm_avg'])[thnq_data==thnq] )
        
        pwiaXsec_JML = np.array( (kin['pwiaXsec_theory'])[thnq_data==thnq] )
        fsiXsec_JML = np.array( (kin['fsiXsec_theory'])[thnq_data==thnq] )

    
        #if(thnq==75):
        #interpolate
        print('pmavg = ',pm_avg)
        print('pwiaXsec = ',pwiaXsec_JML)
        f_pwiaXsec_JML = interp1d(pm_avg, pwiaXsec_JML, fill_value='extrapolate', kind='linear')
        f_fsiXsec_JML = interp1d(pm_avg, fsiXsec_JML, fill_value='extrapolate', kind='linear')
        return pm_avg, f_pwiaXsec_JML, f_fsiXsec_JML
        
        #else:
        #    f_pwiaXsec_JML = interp1d(pm_avg, pwiaXsec_JML, fill_value='extrapolate', kind='cubic')
        #    f_fsiXsec_JML = interp1d(pm_avg, fsiXsec_JML, fill_value='extrapolate', kind='cubic')
        #    return pm_avg, f_pwiaXsec_JML, f_fsiXsec_JML
        
