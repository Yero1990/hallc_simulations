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
import LT.box as B
import matplotlib.pyplot as plt

   
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


def calc_elec_kin(Ei=-1, th_e=-1, kf=-1, Q2=-1, xbj=-1, e_arm='SHMS', case='case1',):


    # Calculates electron kinematics ONLY
    if(case=='case1'):

        # Case 1: Given beam energy (Ei), electron angle (th_e) and momentum (kf)

        # Calculate electron final energy
        Ee = np.sqrt(kf**2 + me**2)

        #4-momentum transfer [GeV^2]
        Q2 = 4. * Ei * Ee * np.sin(th_e/2. * dtr)**2
    
        # Calculate energy transer omega = (Ei - Ee) [GeV]
        omega = Ei - Ee

        # Calculate x-Bjorken
        xbj = Q2 / (2. * MP * omega)
        

        # Calculate electron scattering angle, th_e [deg]
        th_e = 2. * np.arcsin( np.sqrt( Q2 / ( 4 * Ei * Ee ) ) ) / dtr

        print('------------------------------\n'
              'Electron Kinematics, given:  \n'
              'Beam Energy (Eb)          = %.3f GeV' % (Ei),'\n'
              'e- angle (th_e)           = %.3f deg' % (th_e),'\n'            
              'e- momentum (kf)          = %.3f GeV/c' % (kf),'\n'
              '------------------------------\n'
              '4-momentum transfer (Q2)  = %.3f (GeV/c)^2' % (Q2),'\n'
              'x-Bjorken                 = %.3f'    % (xbj),'\n'
              'e- final energy (Ee)      = %.3f GeV' % (Ee),'\n'
              'energy transfer           = %.3f GeV' % (omega),'\n'     
              )
    
    # Calculates electron kinematics ONLY
    if(case=='case4'):

        # Case 4: Given Q2, xbj, and beam energy (Ei) 

        # Calculate energy transer omega = (Ei - Ef) [GeV]
        omega = Q2 / (2. * MP * xbj)    

        # Calculate final electron energy and momentum [GeV]
        Ee = Ei - omega
        kf = np.sqrt(Ee**2 - me**2)

        # Calculate electron scattering angle, th_e [deg]
        th_e = 2. * np.arcsin( np.sqrt( Q2 / ( 4 * Ei * Ee ) ) ) / dtr

        print('------------------------------\n'
              'Electron Kinematics, given:  \n'
              'Beam Energy (Eb)          = %.3f GeV' % (Ei),'\n'
              '4-momentum transfer (Q2)  = %.3f (GeV/c)^2' % (Q2),'\n'
              'x-Bjorken                 = %.3f'    % (xbj),'\n'
              '------------------------------\n'
              'e- angle (th_e)           = %.3f deg' % (th_e),'\n'            
              'e- momentum (kf)          = %.3f GeV/c' % (kf),'\n'
              'e- final energy (Ee)      = %.3f GeV' % (Ee),'\n'
              'energy transfer           = %.3f GeV' % (omega),'\n'     
              )
        
        
# Main Functions (Call other functions above)
if __name__ == "__main__":

    print('Main')

    # Call function to calculate H(e,e'p) elastic kinematics
    #kf, th_p, Pf, Q2, q, omega, xbj = calc_h2(10.6, 12.1686, -1, 'SHMS', 'case1', True)

    
    calc_elec_kin(Ei=10.6, th_e=-1, kf=-1, Q2=2.1, xbj=0.976, e_arm='SHMS', case='case4')
