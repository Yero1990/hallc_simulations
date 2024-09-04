'''
Script to calculate the deuteron d(e,e'p)n reaction kinematics
to determine at which central settings (in angle and momentum)
the Hall C spectrometers should be placed

Author: C. Yero
email: cyero@jlab.org
'''

import numpy as np
#import LT.box as B
import matplotlib.pyplot as plt
#import pandas as pd




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
    Ei = 8.8  #beam energy

    #Q2 = 2.9   #4-momentum transfer ( this can be ignorde for now)

    #Set Q2 Range to cover [GeV^2]
    Q2_min = 1.49 #2.9
    Q2_step = 0.1 
    Q2_max = 1.51 #4.5 + Q2_step   #include endpoint (+Pr_step)
    Q2_range = np.arange(Q2_min, Q2_max, Q2_step)
    
    #Set Missing Momentum Range to cover [GeV]
    Pr_min = 0.1
    Pr_step = 0.05    
    Pr_max = 0.5 + Pr_step   #include endpoint (+Pr_step)
    Pr_range = np.arange(Pr_min, Pr_max, Pr_step)
    
    #Set x-Bjorken Range to cover
    xbj_min = 1.0
    xbj_step = 0.05
    xbj_max = 2. + xbj_step
    xbj_range = np.arange(xbj_min, xbj_max, xbj_step)
    
    #output file to write kinematics
    #fname = 'polarized_deut_kin_summary_Eb%.2f_phi180.csv' % (Ei)
    #fname = 'polarized_deut_kin_summary_Eb%.2f_phi180_HMSwideOpen_thrq35.txt' % (Ei)
    fname = 'd2pol_Eb8p8_phi0.txt'
    
    ofile = open(fname, 'w')
    ofile.write('# d(e,e\'p)n Central Kinematics Summary\n')
    ofile.write('# Beam Energy (Ei) = %.3f GeV\n' % (Ei))
    ofile.write('# 4-Momentum Transfer (Q2) = %.2f - %.2f GeV (step: %.2f) \n' % (Q2_min, Q2_max-Q2_step, Q2_step))
    ofile.write('# x-Bjorken (xbj) = %.2f - %.2f (step: %.2f) \n' % (xbj_min, xbj_max, xbj_step))    
    ofile.write('# Missing Momentum (Pr) = %.2f - %.2f GeV (step: %.2f) \n' % (Pr_min, Pr_max, Pr_step))
    ofile.write('# Hadron Out-of-Plane Angle (phi): rotation axis is q-vector \n')
    ofile.write('# thp = thq - thpq, phi = 0  (q-vector scatters at smaller  angles than proton scattering angle)')
    ofile.write('# \n')
    ofile.write('# ')
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
    
    #ofile.write('#! Pr[f,0]/ \t  xbj[f,1]/ \t kf[f,2]/ \t th_e[f,3]/ \t Pf[f,4]/ \t th_p[f,5]/ \t q[f,6]/ \t th_q[f,7]/ \t th_nq[f,8]/ \t th_pq[f,9]/ \t Q2[f,10]/\n')
    ofile.write('Pr,xbj,kf,th_e,Pf,th_p,q,th_q,th_nq,th_pq,Q2\n')


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
                #thp = thq + thpq  # phi = 180  (q-vector scatters at smaller angles than proton scattering angle)
                thp = thq - thpq   # phi = 0   (q-vector scatters at larger  angles than proton scattering angle)

                if (np.isnan(thp)): continue

                # restrict the proton angle to < 35 deg (allowed by magnet used in polarization)
                #if(thp>=57.5): continue
                

                # restrict the neutron recoil angle relative to q-vector, theta_rq
                #if(thnq < 10. or thnq > 30): continue
                
                if(th_e<7.5): continue
                if(th_e>=40.): continue
                
                ofile.write("  %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f\n" % (Pr, xbj, kf, th_e, Pf, thp, q, thq, thnq, thpq, Q2 ) )
                #ofile.write("%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n" % (Pr, xbj, kf, th_e, Pf, thp, q, thq, thnq, thpq, Q2 ) )
            
    ofile.close()



'''
#plot kinematic correlations
def plot_kin():

    csv_file='polarized_deut_kin_summary_Eb10.55_phi180.csv'
    df           = pd.read_csv(csv_file, comment='#')
    df.to_numpy()

    clr  = ['grey', 'c',    'm',   'r',   'g',   'b', 'darkorange', 'violet', 'gold', 'lightcoral'] #, 'olive', 'sandybrown'] #'darkgray']

    Pm_c = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    fig, ax  = plt.subplots(3, 3)    

    for i in np.arange(len(Pm_c)):    

        print(i)
        print('color=',clr[i])
        print('Pmc=',Pm_c[i])
        
        ax[0,0].plot( df['th_e'].loc[df.Pr == Pm_c[i]] , df['kf'].loc[df.Pr == Pm_c[i]], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
        ax[0,0].set(xlabel='SHMS Angle [deg]', ylabel='SHMS Momentum [GeV/c]')

        
        # Q2 vs. th_e
        ax[0,1].plot( df['th_e'] , df['Q2'], linestyle='', marker='o', mec='k', mfc=clr[1], alpha=0.8)
        ax[0,1].set(xlabel='SHMS Angle [deg]', ylabel=r'$Q^{2}$ [GeV$^{2}$]')
        
        # xbj vs. th_e
        ax[0,2].plot( df['th_e'] , df['xbj'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
        ax[0,2].set(xlabel='SHMS Angle [deg]', ylabel=r'$x_{Bj}$')
        
        # Pr vs. th_e
        ax[1,0].plot( df['th_e'] , df['Pr'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
        ax[1,0].set(xlabel='SHMS Angle [deg]', ylabel=r'Recoil Momentum, $P_{r}$ [GeV/c]')
        
        # th_nq vs. th_e
        #ax[1,1].plot( df['th_e'] , df['th_nq'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
        #ax[1,1].set(xlabel='SHMS Angle [deg]', ylabel=r'$\theta_{rq}$ [deg]')
        
        # th_nq vs. th_e
        ax[1,1].plot( df['th_p'] , df['th_nq'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
        ax[1,1].set(xlabel='HMS Angle [deg]', ylabel=r'$\theta_{rq}$ [deg]')
        
        # Pf vs th_e
        ax[1,2].plot( df['th_e'], df['Pf'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
        ax[1,2].set(xlabel='SHMS Angle [deg]', ylabel='HMS Momentum [GeV/c]')
        
        # Pf vs th_p
        ax[2,0].plot( df['th_p'], df['Pf'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
        ax[2,0].set(xlabel='HMS Angle [deg]', ylabel='HMS Momentum [GeV/c]')
        
        # Q2 vs xbj
        ax[2,1].plot( df['xbj'], df['Q2'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
        ax[2,1].set(xlabel=r'$x_{bj}$', ylabel=r'$Q^{2}$ [GeV$^{2}$]')
        
        # Pr vs th_nq
        ax[2,2].plot( df['th_nq'], df['Pr'], linestyle='', marker='o', mec='k', mfc=clr[i], alpha=0.8)
        ax[2,2].set(xlabel=r'$\theta_{rq}$ [deg]', ylabel=r'Recoil Momentum, $P_{r}$ [GeV/c]')
        
        
        
    plt.tight_layout()
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.show()


'''

if __name__ == "__main__":
    calc_d2_kin()

    #plot_kin()

