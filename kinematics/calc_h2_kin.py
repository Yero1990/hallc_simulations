'''
Script to calculate the hydrogen H(e,e'p) elastics reaction kinematics
to determine at which central settings (in angle and momentum)
the Hall C spectrometers should be placed

Author: C. Yero
email: cyero@jlab.org
'''

import numpy as np

#degrees-to-radians
dtr = np.pi / 180.
    
#Define Particle Mass [GeV]
MP = 0.938272;  #proton
MN = 0.939565;  #neutron
MD = 1.875612;  #deuteron
me = 0.000511;  #electron

    
import numpy as np
#import LT.box as B
#import matplotlib.pyplot as plt




def calc_h2_kin():
    print('main')

  
    #For H(e,e'p) Elastics, require Xbj = 1
    #Recall: Xbj = Q2 / ( 2Mp(E-E') )
    

    #Read deuteron kinematics settings to be used in calculating the corresponding H(e,e'p) elastics 
    fname = 'd2kin_summary_Eb10.60_reduced.txt'
    f = B.get_file(fname)
    Pr = B.get_data(f, 'Pr')      # recoil momentum setting in deuteron [GeV]
    the = B.get_data(f, 'th_e')  #e- scattering angle [deg]


    Ei = 10.6;  #beam energy [GeV]

    #output file to write kinematics
    foname = 'h2_kin_summary_Eb%.2f.txt' % (Ei)
    ofile = open(foname, 'w')
    ofile.write('# H(e,e\'p) Elastics Central Kinematics Summary\n')
    ofile.write('# Beam Energy (Ei) = %.3f GeV\n' % (Ei))
    ofile.write('# \n')
    ofile.write('#! Pr[f,0]/ \t  xbj[f,1]/ \t kf[f,2]/ \t th_e[f,3]/ \t Pf[f,4]/ \t th_p[f,5]/ \t q[f,6]/ \t Q2[f,7]/\n')
    
    #------------------------------------------------------------
    # Case 1: For a fixed beam energy (E), and e- angle (th_e),
    # What does final e- energy (Ee) need to be for Xbj = 1 ?
    #------------------------------------------------------------


    for idx, th_e in enumerate(the):
        #final e- energy required for elastics (given fixed th_e and Ei)
        Ee = 2. * MP * Ei / ( 4. * Ei * np.sin(th_e/2. * dtr)**2 + 2.*MP )  

        #Calculate initial/final e- momentum [GeV / c]
        ki = np.sqrt(Ei**2 - me**2)  #basically ki = Ei, since me is 1/2 MeV (but done, for completeness)
        kf = np.sqrt(Ee**2 - me**2)
        
        #energy transfer [GeV]
        omega = Ei -  Ee 

        #4-momentum transfer [GeV^2]
        Q2 = 4. * Ei * Ee * np.sin(th_e/2. * dtr)**2
        
        xbj = Q2 / (2.* MP * omega)

        print('Ee = ',Ee,' omega = ',omega,' th_e =', th_e, 'Q2 = ', Q2)

        #Calculate 3-momentum transfer magnitude (|q|)
        q = np.sqrt(Q2 + omega**2)
        
        #Calculate final (elastic) proton energy / momentum
        Ef = omega + MP
        Pf = np.sqrt(Ef**2 - MP**2)
        
        
        
        #Calculate final (elastic) proton scattering angle
        sthp = kf * np.sin(th_e * dtr) / Pf  #sine of proton scattering angle
        th_p = np.arcsin(sthp) / dtr  #proton scattering angle [deg] 
        

        #Write to file
        ofile.write("  %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f\n" % (Pr[idx], xbj, kf, th_e, Pf, th_p, q, Q2 ) )
            
    ofile.close()
   
def calc_h2(Ei=-1, th_e=-1, kf=-1, e_arm='SHMS', case='case1', verbose=False):

    if(case=='case1'):
        
        #------------------------------------------------------------
        # Case 1: For a fixed beam energy (E), and e- angle (th_e),
        # What does final e- energy (Ee) need to be for Xbj = 1 ?
        # User Input: let kf = -1, and put actual Ei and th_e
        #------------------------------------------------------------

        #final e- energy required for elastics (given fixed th_e and Ei)
        Ee = 2. * MP * Ei / ( 4. * Ei * np.sin(th_e/2. * dtr)**2 + 2.*MP )  

        #Calculate initial/final e- momentum [GeV / c]
        ki = np.sqrt(Ei**2 - me**2)  #basically ki = Ei, since me is 1/2 MeV (but done, for completeness)
        kf = np.sqrt(Ee**2 - me**2)

    if(case=='case2'):
        
        #----------------------------------------------------------------
        # Case 2: For a fixed beam energy (E), and e- final momentum (kf),
        # What does final e- angle (th_e) need to be for Xbj = 1 ?
        # User Input: let th_e = -1, and put actual Ei and kf
        #----------------------------------------------------------------

        #Calculate initial/final e- momentum [GeV / c]
        ki = np.sqrt(Ei**2 - me**2)  #basically ki = Ei, since me is 1/2 MeV (but done, for completeness)
        Ee = np.sqrt(kf**2 + me**2)

        #calculate e- angle [deg] (given fixed Ei and kf)
        th_e = 2. * np.arcsin( np.sqrt( 1. / (4.*Ei) * (2.* MP * Ei / Ee - (2 * MP) ) ) ) / dtr

    
    #energy transfer [GeV]
    omega = Ei -  Ee 
    
    #4-momentum transfer [GeV^2]
    Q2 = 4. * Ei * Ee * np.sin(th_e/2. * dtr)**2
    
    xbj = Q2 / (2.* MP * omega)
    
    print('Ee = ',Ee,' omega = ',omega,' th_e =', th_e, 'Q2 = ', Q2)
    
    #Calculate 3-momentum transfer magnitude (|q|) [GeV/c]
    q = np.sqrt(Q2 + omega**2)
    
    #Calculate final (elastic) proton energy / momentum
    Ef = omega + MP
    Pf = np.sqrt(Ef**2 - MP**2)
            
    #Calculate final (elastic) proton scattering angle
    sthp = kf * np.sin(th_e * dtr) / Pf  #sine of proton scattering angle
    th_p = np.arcsin(sthp) / dtr  #proton scattering angle [deg]
    
    
    # Spectrometer Momentum Coverage Calcularions
    if(e_arm=='SHMS'):
        
        #e-arm momentum acceptance delta = (P - Pcentral)/Pcentral * 100 [%] 
        dmin_e = -10  
        dmax_e = 22
        
        #hadron arm (HMS)
        dmin_h = -9
        dmax_h = 9
        h_arm = 'HMS'
        
    #calculate minimum/maximum possible particle momentum P, given the spec. central momentum 
    pmin_e = kf * (1. + dmin_e/100.)
    pmax_e = kf * (1. + dmax_e/100.)

    pmin_h = Pf * (1 + dmin_h/100.) 
    pmax_h = Pf * (1 + dmax_h/100.) 
        
    if(verbose):
        print('------------------------------\n'
              'H(e,e\'p) Elastics for  \n'
              'Beam Energy (Eb) = %.3f GeV' % (Ei),'\n'
              'e- angle (th_e)  = %.3f deg' % (th_e),'\n'
              '------------------------------\n'
              'e- momentum (kf)          = %.3f GeV/c' % (kf),'\n'
              'proton angle (th_p)       = %.3f deg' % (th_p),'\n'
              'proton momentum (Pf)      = %.3f GeV/c' % (Pf),'\n'
              '4-momentum transfer (Q2)  = %.3f (GeV/c)^2' % (Q2),'\n'
              '3-momentum transfer (|q|) = %.3f GeV/c' % (q),'\n'
              'energy transfer           = %.3f GeV' % (omega),'\n'
              'x-Bjorken                 = %.3f'    % (xbj),'\n'
              '\n'
              '------------------------------\n'
              '%s (e- arm) Momentum Acceptance        ' % (e_arm),'\n'
              '------------------------------\n'
              '(delta_min, delta_max) = %.1f, %.1f ' % (dmin_e, dmax_e),'\n'
              '(pmin, pmax) = %.3f, %.3f GeV/c ' % (pmin_e, pmax_e),'\n'            
              '------------------------------\n'
              '%s (h arm) Momentum Acceptance        ' % (h_arm),'\n'
              '------------------------------\n'
              '(delta_min, delta_max) = %.1f, %.1f ' % (dmin_h, dmax_h),'\n'
              '(pmin, pmax) = %.3f, %.3f GeV/c ' % (pmin_h, pmax_h),'\n'       
              )
              
    return (kf, th_p, Pf, Q2, q, omega, xbj)  

    
if __name__ == "__main__":

    print('Main')

    
    #calc_h2_kin()

    kf, th_p, Pf, Q2, q, omega, xbj = calc_h2(10.6, 12.1686, -1, 'SHMS', 'case1', True)

