'''
Script to calculate the deuteron d(e,e'p)n reaction kinematics
to determine at which central settings (in angle and momentum)
the Hall C spectrometers should be placed

Author: C. Yero
email: cyero@jlab.org
'''

import numpy as np
import LT.box as B
import matplotlib.pyplot as plt




def calc_d2_kin():
    print('main')

    #degrees-to-radians
    dtr = np.pi / 180.
    
    #Define Particle Mass [GeV]
    MP = 0.938272  #proton
    MN = 0.939565  #neutron
    MD = 1.875612  #deuteron
    MC12 = 11.1878988
    me = 0.000511  #electron

    #Initial parameter kinematics [GeV]
    Ei = 10.549  #beam energy

    Q2 = 3.5   #4-momentum transfer ( this can be ignorde for now)

    #Set Q2 Range to cover [GeV^2]
    Q2_min = 2.1
    Q2_step = 0.25    
    Q2_max = 2.1 + Q2_step   #include endpoint (+Pr_step)
    Q2_range = np.arange(Q2_min, Q2_max, Q2_step)
    
    #Set Missing Momentum Range to cover [GeV]
    Pr_min = 0.
    Pr_step = 0.1    
    Pr_max = 1.0 + Pr_step   #include endpoint (+Pr_step)
    Pr_range = np.arange(Pr_min, Pr_max, Pr_step)
    
    #Set x-Bjorken Range to cover
    xbj_min = 1.3
    xbj_step = 0.01
    xbj_max = 1.6 + xbj_step
    xbj_range = np.arange(xbj_min, xbj_max, xbj_step)
    
    #output file to write kinematics
    #fname = 'polarized_deut_kin_summary_Eb%.2f_xbj_%.2f_Q2_%.2f.txt' % (Ei, xbj_min, Q2_min)
    fname = 'polarized_deut_kin_summary_Eb%.2f_Q2_2.1_thrq35.txt' % (Ei)

    ofile = open(fname, 'w')
    ofile.write('# d(e,e\'p)n Central Kinematics Summary\n')
    ofile.write('# Beam Energy (Ei) = %.3f GeV\n' % (Ei))
    ofile.write('# 4-Momentum Transfer (Q2) = %.2f - %.2f GeV (step: %.2f) \n' % (Q2_min, Q2_max-Q2_step, Q2_step))
    ofile.write('# x-Bjorken (xbj) = %.2f - %.2f (step: %.2f) \n' % (xbj_min, xbj_max, xbj_step))    
    ofile.write('# Missing Momentum (Pr) = %.2f - %.2f GeV (step: %.2f) \n' % (Pr_min, Pr_max, Pr_step))
    ofile.write('# \n')
    ofile.write('# th_nq = 35 +/- 1 deg')
    ofile.write('# \n'
                '# Tensor Polarized Deuterium: d(e,e\'p) kinematics :\n'
                '# \n'
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


    #Loop over Neutron Recoil ("Missing") Momentum
    for Pr in Pr_range:
        
        #Loop over 4-Momentum Transfer Q^2
        for Q2 in Q2_range:
    
            #Loop over x-Bjorken Scale
            for xbj in xbj_range:


                #Calculate energy transfer
                omega = Q2 / (2. * xbj * MP)
                
                #Calculate final e- energy / momentum
                Ee = Ei - omega
                kf = np.sqrt(Ee**2 - me**2)
                #Calculate initial e- momentum
                ki = np.sqrt(Ei**2 - me**2)
                
                #Calculate e- scattering angle
                th_e = 2. * np.arcsin( np.sqrt( Q2/(4.*Ee*Ei) ) ) / dtr
                
                #Calculate 3-momentum transfer magnitude (|q|)
                q = np.sqrt(Q2 + omega**2)
                
                #Calculate Final Neutron ("Recoil") Energy
                Er = np.sqrt(MN**2 + Pr**2)
                #Calculate Final Proton Energy
                Ef = omega + MD - Er
                                
                #Calculate Final Proton Momentum
                Pf = np.sqrt(Ef**2 - MP**2)
                
                
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

                #if(thp>50): continue
                #if(thp>35): continue
                if(thnq<34. or thnq>36.): continue
                #if(Q2>Q2_min): continue
                #if(xbj>xbj_min): continue
                #if(Pr!=0.120): continue
                ofile.write("  %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f\n" % (Pr, xbj, kf, th_e, Pf, thp, q, thq, thnq, thpq, Q2 ) )
            
    ofile.close()
            
if __name__ == "__main__":
    calc_d2_kin()

