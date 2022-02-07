'''
Script to calculate A(e,e'p)n reaction kinematics
to determine at which central settings (in angle and momentum)
the Hall C spectrometers should be placed

Author: C. Yero
email: cyero@jlab.org
'''

import numpy as np
import LT.box as B
import matplotlib.pyplot as plt




def calc_eepX_kin():
    print('main')

    #degrees-to-radians
    dtr = np.pi / 180.

    # 1 atomic mass unit conversion to GeV/c^2
    amu = 0.93149432 # GeV/c^2

    # CaFe potential targets
    # d, 4He, 9Be, 10B, 11B, 12C, 28Si, 40Ar, 48Ti, 40Ca, 48Ca, and 54Fe 

    # Source of atomic mass : https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    particle_dict = {
        'e'   : 5.48579909070e-4 * amu,  # electron
        'n'   : 1.00866491588 * amu,     # neutron  
        '1H'  : 1.00727646658 * amu,
        '2H'  : 2.01410177812 * amu,
        '4He' : 4.002603254   * amu,
        '9Be' : 9.012183065   * amu,
        '10B' : 10.01293695   * amu,
        '11B' : 11.00930536   * amu,
        '12C' : 12.00000      * amu,
        '28Si': 27.976926534  * amu,
        '40Ar': 39.9623831237 * amu,
        '48Ti': 47.94794198   * amu,
        '40Ca': 39.962590863  * amu,
        '48Ca': 47.95252276   * amu,
        '54Fe': 53.93960899   * amu
    }

    #Define Particle Mass [GeV]
    MP = 0.938272  #proton
    MN = 0.939565  #neutron
    MD = 1.875612  #deuteron
    me = 0.000511  #electron
    


    # General masses for Target(e,e'p)X  - one-proton knock-out
    me = particle_dict['e']    # e-beam beam particle mass
    Mt = particle_dict['12C']   # target mass
    Mx = particle_dict['1H']   # detected particle mass
    Mr = particle_dict['11B']    # recoil system mass

    #Initial parameter kinematics [GeV]
    Ei = 10.6  #beam energy
    
    #output file to write kinematics
    fname = 'CaFe_kin_summary_Eb%.2f_C12_MF_kin1.txt' % (Ei)
    ofile = open(fname, 'w')
    ofile.write('# (e,e\'p)X Central Kinematics Summary\n')
    ofile.write('# Beam Energy (Ei) = %.3f GeV\n' % (Ei))
    ofile.write('# \n'
                '# Header Definitions: \n'
                '# Pr   : Central Missing Momentum [GeV/c] \n'
                '# xbj  : Bjorken x  \n'
                '# kf   : e- final momentum [GeV/c] \n'
                '# th_e : e- scattering angle [deg] \n'
                '# Pf   : final proton momentum [GeV/c] \n'
                '# th_p : proton scattering angle [deg] \n'
                '# q    : 3-momentum transfer [GeV] \n'
                '# th_q : relative angle between q-vector and +z (lab)  [deg] \n'
                '# th_nq: relative angle between q-vector and recoil neutron [deg] \n'
                '# th_pq: relative angle between q-vector and knocked-out proton [deg] \n'
                '# Q2   : 4-Momentum Transfer [GeV^2]  \n'
                '# \n'
                )
    
    ofile.write('#! Pr[f,0]/ \t  xbj[f,1]/ \t kf[f,2]/ \t th_e[f,3]/ \t Pf[f,4]/ \t th_p[f,5]/ \t q[f,6]/ \t th_q[f,7]/ \t th_nq[f,8]/ \t th_pq[f,9]/ \t Q2[f,10]/\n')


    # ===================================
    #  Dien' Notes (Starting values):

    # Q2 = 2.1
    # x = 1.2
    # th_e = 8.2 deg
    # pE = 9.67 GeV/c
    # Pr: 140 - 500 MeV/c


    # 02/07/22 | CaFe Open Kinematics to start optimization (to determine which config has the highest rates)
    # kin1:  Eb = 10.6 GeV, Theta_e = 7.5 degree, Pe = 8.55 GeV/c (Q2 = 1.55 GeV^2,   x = Q2/(2*Mp*nu) ~ 0.4 
    # Kin2:  Eb = 10.6 GeV, Theta_e = 6.5 degree, Pe = 8.55 GeV/c (Q2 = 1.16 GeV^2,   x = Q2/(2*Mp*nu) ~ 0.3 
    # kin3:  Eb = 8.6 GeV,  Theta_e = 9.5 degree, Pe = 8 GeV/c (Q2 = 1.88 GeV^2,      x = Q2/(2*Mp*nu) ~ 1.6
    
    
    # ===================================
    
    #Set Q2 Range to cover [GeV^2]
    Q2_min = 0.1
    Q2_step = 0.01    
    Q2_max = 2 + Q2_step   #include endpoint (+Pr_step)
    Q2_range = np.arange(Q2_min, Q2_max, Q2_step)
    
    #Set Missing Momentum Range to cover [GeV]
    Pr_min = 0.150
    Pr_step = 0.01    
    Pr_max = 0.150 + Pr_step   #include endpoint (+Pr_step)
    Pr_range = np.arange(Pr_min, Pr_max, Pr_step)
    
    #Set x-Bjorken Range to cover
    xbj_min = 0.1
    xbj_step = 0.001
    xbj_max = 1 + xbj_step
    xbj_range = np.arange(xbj_min, xbj_max, xbj_step)
    
    #Loop over 4-Momentum Transfer Q^2
    for Q2 in Q2_range:
    
        #Loop over Recoil ("Missing") Momentum
        for Pr in Pr_range:

            #Loop over x-Bjorken Scale
            for xbj in xbj_range:


                #Calculate energy transfer
                omega = Q2 / (2. * xbj * Mx)
                
                #Calculate final e- energy / momentum
                Ee = Ei - omega
                kf = np.sqrt(Ee**2 - me**2)
                #Calculate initial e- momentum
                ki = np.sqrt(Ei**2 - me**2)
                
                #Calculate e- scattering angle
                th_e = 2. * np.arcsin( np.sqrt( Q2/(4.*Ee*Ei) ) ) / dtr
                
                #Calculate 3-momentum transfer magnitude (|q|)
                q = np.sqrt(Q2 + omega**2)
                
                #Calculate Final ("Recoil") Energy
                Er = np.sqrt(Mr**2 + Pr**2)
                #Calculate Final Proton Energy
                Ef = omega + Mt - Er
                #Calculate Final Proton Momentum
                Pf = np.sqrt(Ef**2 - Mx**2)
                
                
                #---------------------------------------
                # Calculate angles relative to q-vector
                #---------------------------------------
                
                #theta_q (relative angle between q-vector of +z (lab) )
                cthq = (ki**2 + q**2 - kf**2) /  (2.*ki*q) #Cosine (theta_q)
                thq = np.arccos(cthq) / dtr   # theta_q [deg]
                
                #theta_pq (relative angle between q-vector and final proton momentum)
                cthpq = (q**2 + Pf**2 - Pr**2) / (2.*q*Pf)
                thpq = np.arccos(cthpq) / dtr  #theta_pq [deg]
          
                #theta_nq (relative angle between q-vector and recoil momentum)
                cthnq = (q**2 + Pr**2 - Pf**2) / (2.*q*Pr)
                thnq = np.arccos(cthnq) / dtr  #theta_nq [deg]
                
                #theta_p (proton angle relative to +z (lab))
                thp = thq + thpq  #this is assuming proton is detected in the forward spec. momentum ( < 90 deg)

                if (np.isnan(thp)): continue
                
                ofile.write("  %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f\n" % (Pr, xbj, kf, th_e, Pf, thp, q, thq, thnq, thpq, Q2 ) )
            
    ofile.close()
    


if __name__ == "__main__":

    calc_eepX_kin()

