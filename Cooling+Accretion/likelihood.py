import numpy as np
import scipy.integrate as integrate
import math
from scipy import interpolate
from detection import *
# from read_SN1987a import*

#tempo_K,tempo_I,tempo_B, E_K,iE_K,theta_K,itheta_K,E_I,iE_I,theta_I,itheta_I,E_B,iE_B = read_SN_csv()
from read_SN1987A_data import tempo_K,tempo_I,tempo_B, E_K,iE_K,theta_K,itheta_K,B_rate_K,E_I,iE_I,theta_I,itheta_I,E_B,iE_B,B_rate_B
  
############### Likelihood ###############
def Likelihood(x, mass_ord, acc):
  # T_c[MeV], T_a[MeV], tau_c[s],tau_a[s], Rc[km], Ma[M_sol]
  Tc0, tau_c, Rc, tau, Ta0, tau_a, Ma, Pee, toff_K, toff_I, toff_B = x
  #print(x)

  Emin = 0.52 #MeV
  Emax = 20*max(Tc0,Ta0)#MeV
  #Emax = 100#MeV
  E_vec=np.linspace(Emin,Emax,100)
  tmin = 0 #s
  tmax = 6*max(tau_c,tau_a) #s
  #tmax = 30 #s
  t_vec=np.linspace(tmin,tmax,100)
  cos_theta=0
  cos_theta_vec = np.linspace(-1, 1, 10)

  #Total number of events  
  ##### Kamiokande-II ##########
  integrand=integrate.simps(dN_Kamiokande(t_vec[:,np.newaxis], E_vec[:,np.newaxis,np.newaxis], cos_theta_vec[:,np.newaxis,np.newaxis,np.newaxis],Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc,eff_typ="total"), cos_theta_vec, axis=0)
  N_tot_Kamiokande= integrate.simps(integrate.simps(integrand, E_vec, axis=0), t_vec, axis=0)[0]
  L_Kamiokande = 2*N_tot_Kamiokande
  
  for i in range(len(theta_K)):
      cos_theta=math.cos(theta_K[i])
      ti = tempo_K[i] + toff_K
      dN=dN_Kamiokande_exp(ti,E_K[i],iE_K[i],cos_theta, Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)[0]
      if dN<=10**(-320):
       dN=10**(-320)
      L_Kamiokande = L_Kamiokande-2*np.log(dN+B_rate_K[i]/2)
  #print("L_K",L_Kamiokande)  
      
  ##### IMB ##########
  integrand=integrate.simps(dN_IMB(t_vec[:,np.newaxis], E_vec[:,np.newaxis,np.newaxis], cos_theta_vec[:,np.newaxis,np.newaxis,np.newaxis],Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc,eff_typ="total"), cos_theta_vec, axis=0)
  N_tot_IMB= integrate.simps(integrate.simps(integrand, E_vec, axis=0), t_vec, axis=0)[0]
  L_IMB = 2*N_tot_IMB
  
  for i in range(len(theta_I)):
      cos_theta=math.cos(theta_I[i])
      ti = tempo_I[i] + toff_I
      dN=dN_IMB_exp(ti,E_I[i],iE_I[i],cos_theta, Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)[0]
      if dN<=10**(-320):
       dN=10**(-320)
      L_IMB = L_IMB-2*np.log(dN)
  #print("L_IMB",L_IMB) 
  
  ##### Baksan ##########
  integrand=integrate.simps(dN_Baksan(t_vec[:,np.newaxis], E_vec[:,np.newaxis,np.newaxis], cos_theta_vec[:,np.newaxis,np.newaxis,np.newaxis],Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc,eff_typ="total"), cos_theta_vec, axis=0)
  N_tot_Baksan= integrate.simps(integrate.simps(integrand, E_vec, axis=0), t_vec, axis=0)[0]
  L_Baksan = 2*N_tot_Baksan
  
  for i in range(len(tempo_B)):
      ti = tempo_B[i] + toff_B
      dN=dN_Baksan_exp(ti,E_B[i],iE_B[i],0, Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)[0]
      if dN<=10**(-320):
       dN=10**(-320)
      L_Baksan = L_Baksan-2*np.log(dN+B_rate_B[i]/2)
  return L_Kamiokande+L_IMB+L_Baksan
Likelihood_vec=np.vectorize(Likelihood)


###############################################################################################################
############################################### Calling Likelihood for NH, IH or Pee ##########################
#  Tc0, tau_c, Rc,tau,Ta0,tau_a, Ma, Pee,mass_ord,acc, toff_K,toff_I,toff_B=x

def Likelihood_c_a_NH(x):
  # Tc0, tau_c, Rc, tau, Ta0, tau_a, Ma, Pee = 0, toff_K, toff_I, toff_B
  x_aux=[x[0], x[1], x[2], x[3], x[4], x[5], x[6], 0, x[7], x[8], x[9]]
  return Likelihood(x_aux,"NH","yes")

def Likelihood_c_a_IH(x):
  # Tc0, tau_c, Rc, tau, Ta0, tau_a, Ma, Pee = 0, toff_K, toff_I, toff_B
  x_aux=[x[0], x[1], x[2], x[3], x[4], x[5], x[6], 0, x[7], x[8], x[9]]
  return Likelihood(x_aux,"IH","yes")

def Likelihood_c_a_Pee(x):
  # Tc0, tau_c, Rc, tau, Ta0, tau_a, Ma, toff_K, toff_I, toff_B, Pee
  x_aux=[x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[10], x[7], x[8], x[9]]
  return Likelihood(x_aux,"no","yes")

def Likelihood_c_a_Pee_cont(x):
  # Tc0, tau_c, Rc, tau, Ta0, tau_a, Ma, toff_K, toff_I, toff_B,  Pee
  x_aux=[x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[10], x[7], x[8],x[9]]
  return Likelihood(x_aux,"no","yes_contemporaneous")

