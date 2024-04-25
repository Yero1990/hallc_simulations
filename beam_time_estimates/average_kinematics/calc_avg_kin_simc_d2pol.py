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



#------------------------------------------------------------
# header information for the output file (FOR JLAB22 GeV / deuteron calculations)
header = \
"""
# Averaged Kinematics Numerical File
#
# NOTE: averaged quantities are determined from SIMC @ the vertex
# calculated quantities are calculated external to SIMC assuming
# (using averged quantities as input, if needed)
#
# id_bin     : 2d bin number
# id_xbin    : xbin number
# id_ybin    : ybin number
# thnq_b     : central theta_nq value (e.g. 5 15 25 35 . . . deg)
# pm_b       : central missing momentum value (80 120 160 . . . MeV)
# Ei         : average beam energy [MeV]
# kf         : average final e- momentum [MeV]
# the        : average final e- angle [deg]
# nu         : average energy transfer [MeV]
# nu_calc    : calculated energy transfer [MeV]
# Q2_calc    : calculated 4-momentum transfer [MeV^2]
# q_calc     : calculated 3-momentum transfer [MeV]
# Ep_calc    : calculated final proton energy [MeV]
# Pf         : average final proton momentum [MeV]
# Pm         : average missing momentum [MeV]
# Pm_calc    : calculated missing momentum [MeV]
# En_calc    : calculated recoil neutron energy [MeV]
# Pf_Par_q   : calculated Pf component (q-frame) parallel to q
# Pf_Perp_q  : calculated Pf component (q-frame) perpendicular to q
# thpq       : avergaed angle between proton and q [deg]
# thpq_calc  : calculated angle between proton and q [deg]
# thpq_cm    : calculated angle between proton and q in c.m. [deg]
# thnq       : averaged angle between recoil neutron and q [deg]
# thnq_calc  : calculated angle between recoil neutron and q [deg]
# cphi_pq    : averaged cos(angle) betweeen proton and q [deg]
# sphi_pq    : averaged sin(angle) betweeen proton and q [deg] 
# alpha_calc : calculated spectatror (neutron) light-cone momentum fraction 
# contz      : bin content for the specified (Pm thnq) bin
# 
id_bin,id_xbin,id_ybin,thnq_b,pm_b,Ei,kf,the,nu,nu_calc,Q2_calc,q_calc,Ep_calc,Pf,Pm,Pm_calc,En_calc,Pf_Par_q,Pf_Perp_q,thpq,thpq_calc,thpq_cm,thnq,thnq_calc,cphi_pq,sphi_pq,alpha_calc,contz
"""
#------------------------------------------------------------





#print argv
#usage: /apps/python/2.7.12/bin/python calc_avg_kin.py 80 fsi 1 Em_final40MeV


# User Set the general file name Eb, Pr, thrq
list_of_args = sys.argv


# deuteron FSI proposal files
#basename='d2_pm800_thrq49_fsi_rad_output' 

# deuteron polarized proposal files
basename='d2_pm400_Q2_2p0_fsi_rad_fieldON_phi0_output'
#basename='c12_pm350_Q2_3p5_fsi_rad_output'
#basename='he4_pm350_Q2_3p5_fsi_rad_output'

output_file = basename+'_avgkin.csv'


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
o.write('# dx = {0:}\n'.format(repr(all.dx)))
o.write('# dy = {0:}\n'.format(repr(all.dy)))
o.write('# nx = {0:}\n'.format(repr(all.nx)))
o.write('# ny = {0:}\n'.format(repr(all.ny)))
o.write('# xmin = {0:}\n'.format(repr(all.xmin)))
o.write('# ymin = {0:}\n'.format(repr(all.ymin)))
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
bin_info_nu        = BI.get_histo_data_arrays(rf.H_nu_2Davg)           # omega, energy transfer
bin_info_xbj       = BI.get_histo_data_arrays(rf.H_xbj_2Davg)          # Xbj, Bjorken
bin_info_Pm        = BI.get_histo_data_arrays(rf.H_Pm_2Davg)           # Missing Momentum
bin_info_thpq      = BI.get_histo_data_arrays(rf.H_theta_pq_2Davg)     # theta_pq [deg]
bin_info_thrq      = BI.get_histo_data_arrays(rf.H_thrq_2Davg)         # theta_nq [deg]
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
   if (acont == 0):
      # skip zero content bins
      continue
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
      
      # calculate electron kinematics from measured, averaged quantities at vertex
      Ef = np.sqrt(kf*kf + me*me)     # final e- energy
      nu_calc = Ei - Ef               # energy transfer

      Q2_calc = 4.*Ei*Ef*np.sin(the*dtr/2.)**2       # calculated 4-momentum transfer
      q_calc = np.sqrt(Q2_calc + nu_calc*nu_calc)    # calculated |q| in the lab frame

      if q_calc==0.:
         # unphysical, skip
         continue
      
      # calculate hadron kinematics
      Ep_calc = np.sqrt( MP**2 + Pf**2)       # calculated final proton energy

      # calculated missing momentum (assuming deuteron mass)
      Pm_calc2 = (nu_calc+MD-Ep_calc)**2 - MN**2
      if (Pm_calc2 < 0.):
         # if unphysical, set it to the averaged Pm obtained from the avg 2D histos
         print ('calculated pm**2 < 0. ', Pm_calc2, ' use Pm_avg : ', Pm)
         Pm_calc = Pm   
      else:
         Pm_calc = np.sqrt ( Pm_calc2 )
         print('Pm_calc = ', Pm_calc)  
      En_calc = np.sqrt(MN**2 + Pm_calc**2);  # calculated recoil neutron energy

      # center of mass motion
      beta_cm = q_calc/(MD+nu_calc)
      gamma_cm = 1./np.sqrt(1. - beta_cm**2)

      # Momentum Components for Proton (in q-frame)
      Pf_Par_q = ( Pf**2 + q_calc**2 - Pm_calc**2)/ (2.*q_calc)
      Pf_Perp_q2 = Pf**2 - Pf_Par_q**2
      if (Pf_Perp_q2 < 0.):
         Pf_Perp_q = Pf*np.sin(dtr*thpq)
         thpq_calc = thpq
      else:
         Pf_Perp_q = np.sqrt(Pf_Perp_q2)
         cthpq = Pf_Par_q/Pf     #Cos(theta_pq)
         thpq_calc = np.arccos(cthpq)/dtr             
      Pf_Par_q_cm = gamma_cm*Pf_Par_q - gamma_cm*beta_cm*Ep_calc   #parallel component of proton in cm

      # proton angle in the cm
      thp_calc_cm = 0.

      if Pf_Par_q_cm == 0. :
         thp_calc_cm = np.pi
      if Pf_Par_q_cm > 0. :
         thp_calc_cm = np.arctan(Pf_Perp_q/Pf_Par_q_cm)
      if Pf_Par_q_cm < 0. :
         thp_calc_cm = np.pi+np.arctan(Pf_Perp_q/Pf_Par_q_cm)

      thpq_cm = thp_calc_cm/dtr

      # calculate angles using calculated Pmiss
      denom = q_calc**2 + Pm_calc**2 - Pf**2
      num = (2.*q_calc*Pm_calc)

      cth_nq = -2.    #Cos(theta_nq)
      thnq_calc = -1.
      if num > 0. : 
         cth_nq = denom/num
         thnq_calc = 0.
      if abs(cth_nq) <=1.:
         thnq_calc = np.arccos(cth_nq)/dtr;
      # calculate alpha
      pz_n = Pm_calc*np.cos(thnq_calc*dtr)
      p_n_minus = En_calc - pz_n
      alpha_calc = p_n_minus/MN



      # for JLab 22 GeV calculation
      l = "%i,%i,%i,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3E\n"%(
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
                                                                                                          # 8 avg. energy transfer
                                                                                                          nu,  \
                                                                                                          # 9 calc.  energy transfer
                                                                                                          nu_calc, \
                                                                                                          # 10 calc. average 4-Momentum transfer
                                                                                                          Q2_calc, \
                                                                                                          # 11 calc. average |q| 3-momentum transfer
                                                                                                          q_calc, \
                                                                                                          # 12 calc. average final proton energy (assume proton mass)
                                                                                                          Ep_calc, \
                                                                                                          # 13  avg. final proton momentum
                                                                                                          Pf, \
                                                                                                          # 14  avg. missing momentum
                                                                                                          Pm, \
                                                                                                          # 15 calc. average Missing momentum  (assume deuteron mass)
                                                                                                          Pm_calc, \
                                                                                                          # 16 calc. final neutron energy 
                                                                                                          En_calc, \
                                                                                                          # 17 calc. Pf component (q-frame) parallel to q
                                                                                                          Pf_Par_q, \
                                                                                                          # 18 calc. Pf component (q-frame) perpendicular to q
                                                                                                          Pf_Perp_q, \
                                                                                                          # 19
                                                                                                          thpq, \
                                                                                                          # 20
                                                                                                          thpq_calc, \
                                                                                                          # 21
                                                                                                          thpq_cm, \
                                                                                                          # 22
                                                                                                          thnq, \
                                                                                                          # 23
                                                                                                          thnq_calc, \
                                                                                                          # 24
                                                                                                          cphi_pq, \
                                                                                                          # 25
                                                                                                          sphi_pq, \
                                                                                                          # 26
                                                                                                          alpha_calc, \
                                                                                                          # 27
                                                                                                          all.cont[i])
       
                                                                         
      o.write(l)
o.close()

