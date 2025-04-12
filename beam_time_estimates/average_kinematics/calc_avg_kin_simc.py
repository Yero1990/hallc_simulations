#!/usr/bin/python


# C.Y. Jan 17, 2023 | This is a copy of the code to calculate the avergaed kinematics, but modified for JLab 22 GeV variables. I've commented out
# some lines or variables not relevant for the calculation, to make the data simpler.

# calculate averaged kinematics from SIMC analysis of 2D Histos (Pm vs. th_nq bins) for JLab 22 GeV

import sys
from sys import argv
import getopt
import os
from LT import datafile
import bin_info2 as BI


#Set proper paths to import ROOT ( I had these paths in ifarm, but am having trouble running the code due to compatibility with ROOT / python)
sys.path.append('./pyroot/')
sys.path.append('/apps/root/PRO/lib/')
sys.path.append('/apps/root/PRO/')

#so far, works with:
#tags/v6-18-04@v6-18-04
# python 3.43
# but still unable to use: from ROOT import * for python3

# to run the code on ifarm:
# 1. source /apps/root/6.18.04/setroot_CUE.csh 
# 2. set proper basename file to be read
# 3. python calc_avg_kin_simc.py  

# do the root operations directly here
import ROOT as R
#from ROOT import *
#from ROOT import *
import numpy as np


# some constants
dtr = np.pi/180.

#MeV
MP = 938.272
MN = 939.566
MD = 1875.6127
me = 0.51099 

'''
#------------------------------------------------------------
# header information for the output file (ORIGINAL)
header = \
"""
# averaged kinematics results
# averaged kinematic varibles used as input to calculate the averaged 
# kinematics: Ei, omega, th_e, pf
#
# variables with _mc attached are from histograms not calculated
# alpha is the spectatror (neutron) alpha
#\\ xb = th_nq
#\\ yb = pm
# current header line:
#! i_b[i,0]/ i_x[i,1]/ i_y[i,2]/ xb[f,3]/ yb[f,4]/ Ei[f,5]/ kf[f,6]/ th_e[f,7]/ omega_mc[f,8]/ omega[f,9]/ Q2[f,10]/ Q2_calc[f,11]/ q_mc[f,12]/ q_lab[f,13]/ Ep_calc[f,14]/ pf[f,15]/ pm_mc[f,16]/ pm[f,17]/ En_calc[f,18]/ beta_cm[f,19]/ gamma_cm[f,20]/ PfPar_q[f,21]/ PfPerp_q[f,22]/ theta_pq[f,23]/ theta_pq_calc[f,24]/ PfPar_cm[f,25]/ th_pq_cm[f,26]/ th_nq_mc[f,27]/ th_nq_calc[f,28]/  cos_phi[f,29]/  sin_phi[f,30]/  alpha_c[f,31]/  GEp[f,32]/   GMp[f,33]/   sigMott[f,34]/   Ksig_cc1[f,35]/  nx[i,36]/ ny[i,37]/ cont[f,38]/        
"""
#------------------------------------------------------------
'''

#------------------------------------------------------------
# header information for the output file (FOR JLAB22 GeV / deuteron calculations)
header = \
"""
# averaged kinematics results
# averaged kinematic varibles used as input to calculate the averaged 
# kinematics: Ei, omega, th_e, pf
#
# variables with _mc attached are from histograms not calculated
# alpha is the spectatror (neutron) alpha
#\\ xb = th_nq
#\\ yb = pm
# current header line:
#! i_b[i,0]/ i_x[i,1]/ i_y[i,2]/ xb[f,3]/ yb[f,4]/ Ei[f,5]/ kf[f,6]/ th_e[f,7]/ omega[f,8]/ Q2_calc[f,9]/ q_lab[f,10]/ Ep_calc[f,11]/ pf[f,12]/ pm_mc[f,13]/ pm[f,14]/ En_calc[f,15]/ beta_cm[f,16]/ gamma_cm[f,17]/ PfPar_q[f,18]/ PfPerp_q[f,19]/ theta_pq[f,20]/ theta_pq_calc[f,21]/ PfPar_cm[f,22]/ th_pq_cm[f,23]/ th_nq_mc[f,24]/ th_nq_calc[f,25]/  cos_phi[f,26]/  sin_phi[f,27]/  alpha_c[f,28]/  nx[i,29]/ ny[i,30]/ cont[f,31]/        
"""
#------------------------------------------------------------





#print argv
#usage: /apps/python/2.7.12/bin/python calc_avg_kin.py 80 fsi 1 Em_final40MeV


# User Set the general file name Eb, Pr, thrq
list_of_args = sys.argv


#Eb=14   #list_of_args[1]
#Pr=1    #list_of_args[2]
#thrq=26 #list_of_args[3]
#basename='d2_Eb%s_Pr%s_thrq%s_norad_output' %(Eb, Pr, thrq)

# deuteron FSI proposal files
basename='d2_pm800_thrq49_fsi_rad_output' 

# deuteron polarized proposal files
#basename='d2_pm350_Q2_3p5_fsi_rad_output'

output_file = basename+'_avgkin.txt'


o = open(output_file,'w')

#Open root file to read avg kin histos
root_file = basename+'.root'

print('root_fil',root_file)
# open ROOTfile
rf = R.TFile(root_file)

# start with 2D yield histo, Fill(Pm, thnq, FullWeight)
all = BI.get_histo_data_arrays(rf.H_Pm_vs_thrq_v) 
# write the necessary header parameters
o.write('# histogram parameters \n')
o.write('#\ dx = {0:}\n'.format(repr(all.dx)))
o.write('#\ dy = {0:}\n'.format(repr(all.dy)))
o.write('#\ nx = {0:}\n'.format(repr(all.nx)))
o.write('#\ ny = {0:}\n'.format(repr(all.ny)))
o.write('#\ xmin = {0:}\n'.format(repr(all.xmin)))
o.write('#\ ymin = {0:}\n'.format(repr(all.ymin)))
# write header
o.write(header)

#Get 2D Histogram Bin Info (Avg. kin)
bin_info_Ei        = BI.get_histo_data_arrays(rf.H_Ein_2Davg)          #inc. beam energy [GeV]
bin_info_kf        = BI.get_histo_data_arrays(rf.H_kf_2Davg)           #final e- momentum [GeV]
bin_info_the       = BI.get_histo_data_arrays(rf.H_the_2Davg)          #final e- angle [deg]
bin_info_Pf        = BI.get_histo_data_arrays(rf.H_Pf_2Davg)           #final p momentum [GeV]
bin_info_thp       = BI.get_histo_data_arrays(rf.H_thp_2Davg)          #final p angle [deg]
bin_info_q         = BI.get_histo_data_arrays(rf.H_q_2Davg)            # |q| momentum transfer
bin_info_thq       = BI.get_histo_data_arrays(rf.H_theta_q_2Davg)      # q-angle with +z beam
bin_info_Q2        = BI.get_histo_data_arrays(rf.H_Q2_2Davg)           # Q2 4-momentum transfer
bin_info_nu        = BI.get_histo_data_arrays(rf.H_nu_2Davg)        # omega, energy transfer
bin_info_xbj       = BI.get_histo_data_arrays(rf.H_xbj_2Davg)          # Xbj, Bjorken
bin_info_Pm         = BI.get_histo_data_arrays(rf.H_Pm_2Davg)           # Missing Momentum
bin_info_thpq      = BI.get_histo_data_arrays(rf.H_theta_pq_2Davg)     # theta_pq [deg]
bin_info_thrq       = BI.get_histo_data_arrays(rf.H_thrq_2Davg)     # theta_nq [deg]
bin_info_cphi_pq   = BI.get_histo_data_arrays(rf.H_cphi_pq_2Davg)      # cos(phi_pq) (-1,1)
bin_info_sphi_pq   = BI.get_histo_data_arrays(rf.H_sphi_pq_2Davg)      # sin(phi_pq) (-1,1)


#Loop over bin number (xbin, ybin)->(th_nq_bin, Pm_bin)
for i,acont in enumerate(all.cont):
   
   # get bin values
   i_bin = all.i[i]
   i_xbin = all.ix[i]
   i_ybin = all.iy[i]
   thnq_b = all.xb[i]
   pm_b = all.yb[i]
   if (acont == -1):
      # skip zero content bins
      #continue
      print('acont = ',acont)
   else:
      
      # convert rad to deg and GeV to MeV 
      Ei        = bin_info_Ei.cont[i]*1000.      
      kf        = bin_info_kf.cont[i]*1000.
      the       = bin_info_the.cont[i]
      Pf        = bin_info_Pf.cont[i]*1000.
      thp       = bin_info_thp.cont[i]
      q         = bin_info_q.cont[i]*1000.
      thq       = bin_info_thq.cont[i]
      Q2        = bin_info_Q2.cont[i]*1.e6
      nu        = bin_info_nu.cont[i]*1000.
      xbj       = bin_info_xbj.cont[i]
      Pm        = bin_info_Pm.cont[i]*1000.
      thpq      = bin_info_thpq.cont[i]
      thnq      = bin_info_thrq.cont[i]
      cphi_pq   = bin_info_cphi_pq.cont[i]
      sphi_pq   = bin_info_sphi_pq.cont[i]

      print("i = ", i)
      print("thnq_b = ",thnq_b)
      print("pm_b = ",  pm_b)
      print("Ei  = ", Ei)
      print("kf = ", kf)
      print("Pm = ", Pm)
      
      # calculate electron kinematics from measured, averaged quantities
      Ef = np.sqrt(kf*kf + me*me)  
      nu_calc = Ei - Ef

      Q2_calc = 4.*Ei*Ef*np.sin(the*dtr/2.)**2
      q_calc = np.sqrt(Q2_calc + nu_calc*nu_calc)    # |q| in the lab frame
      if q_calc==0.:
         #unphysical, skip
         continue
         # calculate hadron kinematics
      Ep = np.sqrt( MP**2 + Pf**2)
      # calculated missing momentum
      Pm_calc2 = (nu_calc+MD-Ep)**2 - MN**2
      if (Pm_calc2 < 0.):
         print ('calculated pm**2 < 0. ', Pm_calc2, ' use Pm_avg : ', Pm)
         Pm_calc = Pm   #set it to the average Pm from 2D histo
      else:
         Pm_calc = np.sqrt ( Pm_calc2 )
         print('Pm_calc = ', Pm_calc)  
      En_calc = np.sqrt(MN**2 + Pm_calc**2);

      # center of mass motion
      beta_cm = q_calc/(MD+nu_calc)
      gamma_cm = 1./np.sqrt(1. - beta_cm**2)

      # Momentum Components for Proton (in q-frame)
      Pf_par = ( Pf**2 + q_calc**2 - Pm_calc**2)/ (2.*q_calc)
      Pf_perp2 = Pf**2 - Pf_par**2
      if (Pf_perp2 < 0.):
         #print 'calculated Pf_perp**2<0. : ', Pf_perp2,' --->   estimate it using theta_pq :', thpq
         #print 'Pf_par = ', Pf_par, ', Pf_perp(pf) = ', Pf*np.sin(dtr*thpq), ', Pf_perp(pm) = ', Pm*np.sin(dtr*thnq)   
         Pf_perp = Pf*np.sin(dtr*thpq)
         th_pq_calc = thpq
      else:
         Pf_perp = np.sqrt(Pf_perp2)
         cthpq = Pf_par/Pf     #Cos(theta_pq)
         th_pq_calc = np.arccos(cthpq)/dtr             
      Pf_par_cm = gamma_cm*Pf_par - gamma_cm*beta_cm*Ep   #parallel component of proton in cm

      # proton angle in the cm
      thp_calc_cm = 0.

      if Pf_par_cm == 0. :
         thp_calc_cm = np.pi
      if Pf_par_cm > 0. :
         thp_calc_cm = np.arctan(Pf_perp/Pf_par_cm)
      if Pf_par_cm < 0. :
         thp_calc_cm = np.pi+np.arctan(Pf_perp/Pf_par_cm)

      theta_pq_cm = thp_calc_cm/dtr

      # calculate angles using calculated Pmiss
      denom = q_calc**2 + Pm_calc**2 - Pf**2
      num = (2.*q_calc*Pm_calc)

      cth_nq = -2.    #Cos(theta_nq)
      theta_nq_calc = -1.
      if num > 0. : 
         cth_nq = denom/num
         theta_nq_calc = 0.
      if abs(cth_nq) <=1.:
         theta_nq_calc = np.arccos(cth_nq)/dtr;
      # calculate alpha
      pz_n = Pm_calc*np.cos(theta_nq_calc*dtr)
      p_n_minus = En_calc - pz_n
      alpha_calc = p_n_minus/MN

      #Calculate the deForest Cross Section Factor  K * sig_cc1,  where K = Pf * Ep , units: ub * MeV^2 / sr^2
      #sig_Mott =  sigMott(kf, the, Q2_calc)   #ub / sr
      #GE_p = GEp(Q2_calc, 'JRA')
      #GM_p = GMp(Q2_calc, 'JRA')
      #Kfact, f_rec, sig_eN, de_Forest = deForest(Ef, Q2_calc, q_calc, Pf, Pm_calc, the, thp, cphi_pq, th_pq_calc, sig_Mott, GE_p, GM_p)
      
      #print('ix=',i_xbin,' iy=',i_ybin,' pm=',Pm_calc,' Kfact=',Kfact,' f_rec=',f_rec,' sig_eN=',sig_eN)
      # write output file

      # original
      '''
      l = "%i %i %i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %i %i %/3E\n"%( \
                                                                                                                                  # 0 2d bin number
                                                                                                                                  i_bin, \
                                                                                                                                  # 1 
                                                                                                                                  i_xbin, \
                                                                                                                                  # 2
                                                                                                                                  i_ybin, \
                                                                                                                                  # 3 central thnq_bin
                                                                                                                                  thnq_b, \
                                                                                                                                  # 4 central pm_bin
                                                                                                                                  pm_b, \
                                                                                                                                  # 5 avg. beam energy 
                                                                                                                                  Ei, \
                                                                                                                                  # 6 avg. e- momentum
                                                                                                                                  kf, \
                                                                                                                                  # 7 avg. e- angle
                                                                                                                                  the, \
                                                                                                                                  # 8 MC average energy transfer
                                                                                                                                  nu, \
                                                                                                                                  # 9 calc. average energy transer
                                                                                                                                  nu_calc, \
                                                                                                                                  # 10 MC average 4-Momentum tranfer
                                                                                                                                  Q2, \
                                                                                                                                  # 11 calc. average 4-Momentum transfer
                                                                                                                                  Q2_calc, \
                                                                                                                                  # 12 MC average |q| 3-momentum tranfer
                                                                                                                                  q, \
                                                                                                                                  # 13 calc. average |q| 3-momentum transfer
                                                                                                                                  q_calc, \
                                                                                                                                  # 14 calc. average final proton energy (assume proton mass)
                                                                                                                                  Ep, \
                                                                                                                                  # 15 MC average final proton momentum
                                                                                                                                  Pf, \
                                                                                                                                  # 16 MC average missing momentum
                                                                                                                                  Pm, \
                                                                                                                                  # 17 calc. average Missing momentum  (assume deuteron mass)
                                                                                                                                  Pm_calc, \
                                                                                                                                  # 18
                                                                                                                                  En_calc, \
                                                                                                                                  # 19
                                                                                                                                  beta_cm, \
                                                                                                                                  # 20
                                                                                                                                  gamma_cm, \
                                                                                                                                  # 21
                                                                                                                                  Pf_par, \
                                                                                                                                  # 22
                                                                                                                                  Pf_perp, \
                                                                                                                                  # 23
                                                                                                                                  thpq, \
                                                                                                                                  # 24
                                                                                                                                  th_pq_calc, \
                                                                                                                                  # 25
                                                                                                                                  Pf_par_cm, \
                                                                                                                                  # 26
                                                                                                                                  theta_pq_cm, \
                                                                                                                                  # 27
                                                                                                                                  thnq, \
                                                                                                                                  # 28
                                                                                                                                  theta_nq_calc, \ 
                                                                                                                                  # 29
                                                                                                                                  cphi_pq, \
                                                                                                                                  # 30
                                                                                                                                  sphi_pq, \
                                                                                                                                  # 31
                                                                                                                                  alpha_calc, \
                                                                                                                                  # 32
                                                                                                                                  GE_p, \
                                                                                                                                  # 33
                                                                                                                                  GM_p, \
                                                                                                                                  # 34
                                                                                                                                  sig_Mott, \
                                                                                                                                  # 35
                                                                                                                                  de_Forest, \
                                                                                                                                  # 36
                                                                                                                                  all.nx, \
                                                                                                                                  # 37
                                                                                                                                  all.ny, \
                                                                                                                                  # 38
                                                                                                                                  all.cont[i])
      
      '''
      '''
      # for JLab 22 GeV calculation / fsi studies
      l = "%i %i %i %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %i %i %.3E\n"%( \
                                                                                                          # 0 2d bin number
                                                                                                          i_bin, \
                                                                                                          # 1 
                                                                                                          i_xbin, \
                                                                                                          # 2
                                                                                                          i_ybin, \
                                                                                                          # 3 central thnq_bin
                                                                                                          thnq_b, \
                                                                                                          # 4 central pm_bin
                                                                                                          pm_b, \
                                                                                                          # 5 avg. beam energy 
                                                                                                          Ei, \
                                                                                                          # 6 avg. e- momentum
                                                                                                          kf, \
                                                                                                          # 7 avg. e- angle
                                                                                                          the,                
                                                                                                          # 8 calc. average energy transer
                                                                                                          nu_calc, \
                                                                                                          # 9 calc. average 4-Momentum transfer
                                                                                                          Q2_calc, \
                                                                                                          # 10 calc. average |q| 3-momentum transfer
                                                                                                          q_calc, \
                                                                                                          # 11 calc. average final proton energy (assume proton mass)
                                                                                                          Ep, \
                                                                                                          # 12 MC average final proton momentum
                                                                                                          Pf, \
                                                                                                          # 13 MC average missing momentum
                                                                                                          Pm, \
                                                                                                          # 14 calc. average Missing momentum  (assume deuteron mass)
                                                                                                          Pm_calc, \
                                                                                                          # 15
                                                                                                          En_calc, \
                                                                                                          # 16
                                                                                                          beta_cm, \
                                                                                                          # 17
                                                                                                          gamma_cm, \
                                                                                                          # 18
                                                                                                          Pf_par, \
                                                                                                          # 19
                                                                                                          Pf_perp, \
                                                                                                          # 20
                                                                                                          thpq, \
                                                                                                          # 21
                                                                                                          th_pq_calc, \
                                                                                                          # 22
                                                                                                          Pf_par_cm, \
                                                                                                          # 23
                                                                                                          theta_pq_cm, \
                                                                                                          # 24
                                                                                                          thnq, \
                                                                                                          # 25
                                                                                                          theta_nq_calc, \
                                                                                                          # 26
                                                                                                          cphi_pq, \
                                                                                                          # 27
                                                                                                          sphi_pq, \
                                                                                                          # 28
                                                                                                          alpha_calc, \                                          
                                                                                                          # 29
                                                                                                          all.nx, \
                                                                                                          # 30
                                                                                                          all.ny, \
                                                                                                          # 31
                                                                                                          all.cont[i])
      '''
            # for JLab 22 GeV calculation / fsi studies
      l = "%i %i %i %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %i %i %.3E\n"%( \
                                                                                                          # 0 2d bin number
                                                                                                          i_bin, 
                                                                                                          # 1 
                                                                                                          i_xbin, 
                                                                                                          # 2
                                                                                                          i_ybin, 
                                                                                                          # 3 central thnq_bin
                                                                                                          thnq_b, 
                                                                                                          # 4 central pm_bin
                                                                                                          pm_b, 
                                                                                                          # 5 avg. beam energy 
                                                                                                          Ei, 
                                                                                                          # 6 avg. e- momentum
                                                                                                          kf, 
                                                                                                          # 7 avg. e- angle
                                                                                                          the,               
                                                                                                          # 8 calc. average energy transer
                                                                                                          nu_calc, 
                                                                                                          # 9 calc. average 4-Momentum transfer
                                                                                                          Q2_calc, 
                                                                                                          # 10 calc. average |q| 3-momentum transfer
                                                                                                          q_calc, 
                                                                                                          # 11 calc. average final proton energy (assume proton mass)
                                                                                                          Ep, 
                                                                                                          # 12 MC average final proton momentum
                                                                                                          Pf, 
                                                                                                          # 13 MC average missing momentum
                                                                                                          Pm, 
                                                                                                          # 14 calc. average Missing momentum  (assume deuteron mass)
                                                                                                          Pm_calc, 
                                                                                                          # 15
                                                                                                          En_calc, 
                                                                                                          # 16
                                                                                                          beta_cm, 
                                                                                                          # 17
                                                                                                          gamma_cm, 
                                                                                                          # 18
                                                                                                          Pf_par, 
                                                                                                          # 19
                                                                                                          Pf_perp, 
                                                                                                          # 20
                                                                                                          thpq, 
                                                                                                          # 21
                                                                                                          th_pq_calc, 
                                                                                                          # 22
                                                                                                          Pf_par_cm, 
                                                                                                          # 23
                                                                                                          theta_pq_cm, 
                                                                                                          # 24
                                                                                                          thnq, 
                                                                                                          # 25
                                                                                                          theta_nq_calc, 
                                                                                                          # 26
                                                                                                          cphi_pq, 
                                                                                                          # 27
                                                                                                          sphi_pq, 
                                                                                                          # 28
                                                                                                          alpha_calc, 
                                                                                                          #29
                                                                                                          all.nx, 
                                                                                                          #30
                                                                                                          all.ny, 
                                                                                                          all.cont[i])

      
       
                                                                         
      o.write(l)
o.close()

