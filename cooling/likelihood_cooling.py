import numpy as np
import scipy.integrate as integrate
import math
from detection import*
from read_SN1987a import*

#tempo_K,tempo_I,tempo_B, E_K,iE_K,theta_K,itheta_K,E_I,iE_I,theta_I,itheta_I,E_B,iE_B = read_SN_csv()
from read_SN1987A_data import tempo_K,tempo_I,tempo_B, E_K,iE_K,theta_K,itheta_K,E_I,iE_I,theta_I,itheta_I,E_B,iE_B
#Maria Laura Costantini, et.al.
B_rate_K=[1.0*10**-5, 5.4*10**-4, 3.1*10**-2, 8.5*10**-3, 5.3*10**-4, 7.1*10**-2, 5.0*10**-6, 1.0*10**-5, 1.0*10**-5, 1.8*10**-2, 4.0*10**-4, 1.4*10**-2, 7.3*10**-2, 5.2*10**-2, 1.8*10**-2, 7.3*10**-2]
#Lamb and Loredo
#B_rate_K=[1.6*10**-5, 1.9*10**-3, 2.9*10**-2, 1.2*10**-2, 2.1*10**-3, 3.7*10**-2, 4.5*10**-5, 8.2*10**-5, 1.5*10**-5, 1.5*10**-2, 1.9*10**-3, 1.6*10**-2, 3.8*10**-2, 2.9*10**-2, 2.8*10**-2, 3.8*10**-2]
#Baksan
B_rate_B=[8.4*10**-4, 1.3*10**-3, 1.2*10**-3, 1.3*10**-3, 1.3*10**-3]

#Gaussian in energy and exp. uncertainty
def Gaussian_E(E, Edata, sigmaData):
  gauss = np.exp(-(Edata - E)**2./2./sigmaData**2)/(sigmaData*math.sqrt(2.*math.pi))
  return gauss
Gaussian_E_vec=np.vectorize(Gaussian_E)

############################################### Time-Dependent ##########################################

def dN_Kamiokande_gaussian_vec(E, t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no") :
  gauss_R = dN_Kamiokande(t, E,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc) * Gaussian_E(E, Epos, delta_Epos)
  return gauss_R
#dN_Kamiokande_gaussian_vec=np.vectorize(dN_Kamiokande_gaussian)

def dN_IMB_gaussian_vec(E, t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no") :
  gauss_R = dN_IMB(t, E,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc) * Gaussian_E(E, Epos, delta_Epos)
  return gauss_R
#dN_IMB_gaussian_vec=np.vectorize(dN_IMB_gaussian)

def dN_Baksan_gaussian_vec(E, t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no") :
  gauss_R = dN_Baksan(t, E,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc) * Gaussian_E(E, Epos, delta_Epos)
  return gauss_R
#dN_Baksan_gaussian_vec=np.vectorize(dN_Baksan_gaussian)

def dN_Kamiokande_exp(t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no") :
	E_max=Epos+5*delta_Epos
	E_min=Epos-5*delta_Epos
	if E_min<0.52:
		E_min=0.52
	E_vec=np.linspace(E_min,E_max,101)
	#dN = integrate.quad(lambda E: dN_Kamiokande_gaussian(E,t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,Ta0,tau_a,Ma,Pee,mass_ord,acc),0.52, Epos+5*delta_Epos)
	dN = integrate.simps(dN_Kamiokande_gaussian_vec(E_vec[:,np.newaxis], t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc), E_vec, axis=0)
	return dN
	
def dN_IMB_exp(t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no") :
	E_max=Epos+5*delta_Epos
	E_min=Epos-5*delta_Epos
	if E_min<0.52:
		E_min=0.52
	E_vec=np.linspace(E_min,E_max,101)
	dN = integrate.simps(dN_IMB_gaussian_vec(E_vec[:,np.newaxis], t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc), E_vec, axis=0)
	return dN
	
def dN_Baksan_exp(t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no") :
	E_max=Epos+5*delta_Epos
	E_min=Epos-5*delta_Epos
	if E_min<0.52:
		E_min=0.52
	E_vec=np.linspace(E_min,E_max,101)
	dN = integrate.simps(dN_Baksan_gaussian_vec(E_vec[:,np.newaxis], t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc), E_vec, axis=0)
	return dN
dN_Baksan_exp_vec=np.vectorize(dN_Baksan_exp)
  
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
  integrand=integrate.simps(dN_Kamiokande(t_vec[:,np.newaxis], E_vec[:,np.newaxis,np.newaxis], cos_theta_vec[:,np.newaxis,np.newaxis,np.newaxis],Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc), cos_theta_vec, axis=0)
  #integrand=2*dN_Kamiokande_vec(t_vec[:,np.newaxis], E_vec[:,np.newaxis,np.newaxis],0,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)
  N_tot_Kamiokande= integrate.simps(integrate.simps(integrand, E_vec, axis=0), t_vec, axis=0)[0]
  #print("N_K",N_tot_Kamiokande)
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
  integrand=integrate.simps(dN_IMB(t_vec[:,np.newaxis], E_vec[:,np.newaxis,np.newaxis], cos_theta_vec[:,np.newaxis,np.newaxis,np.newaxis],Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc), cos_theta_vec, axis=0)
  #integrand=2*dN_IMB_vec(t_vec[:,np.newaxis], E_vec[:,np.newaxis,np.newaxis],0,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)
  N_tot_IMB= integrate.simps(integrate.simps(integrand, E_vec, axis=0), t_vec, axis=0)[0]
  #print("N_IMB",N_tot_IMB)
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
  integrand=integrate.simps(dN_Baksan(t_vec[:,np.newaxis], E_vec[:,np.newaxis,np.newaxis], cos_theta_vec[:,np.newaxis,np.newaxis,np.newaxis],Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc), cos_theta_vec, axis=0)
  #integrand=2*dN_Baksan_vec(t_vec[:,np.newaxis], E_vec[:,np.newaxis,np.newaxis],0,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)
  N_tot_Baksan= integrate.simps(integrate.simps(integrand, E_vec, axis=0), t_vec, axis=0)[0]
  #print("N_Baksan",N_tot_Baksan)
  L_Baksan = 2*N_tot_Baksan
  
  for i in range(len(tempo_B)):
      #cos_theta=1
      #dN=integrate.simps(dN_Baksan_exp_vec(tempo_B[i],E_B[i],iE_B[i],cos_theta_vec[:,np.newaxis], Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc), cos_theta_vec, axis=0)[0]
      ti = tempo_B[i] + toff_B
      dN=dN_Baksan_exp(ti,E_B[i],iE_B[i],0, Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)[0]
      if dN<=10**(-320):
       dN=10**(-320)
      L_Baksan = L_Baksan-2*np.log(dN+B_rate_B[i]/2)
  #print("L_Baksan",L_Baksan)       
  return L_Kamiokande+L_IMB+L_Baksan
Likelihood_vec=np.vectorize(Likelihood)


def Likelihood_c_NH(x):
  # Tc0, tau_c, Rc, tau, Ta0=0, tau_a=0, Ma=0, Pee=0, toff_K, toff_I, toff_B
  x_aux=[x[0], x[1], x[2], x[3], 0, 0, 0, 0, x[4], x[5], x[6]]
  return Likelihood(x_aux,"NH","no")

def Likelihood_c_IH(x):
  # Tc0, tau_c, Rc, tau, Ta0=0, tau_a=0, Ma=0, Pee=0, toff_K, toff_I, toff_B
  x_aux=[x[0], x[1], x[2], x[3], 0, 0, 0, 0, x[4], x[5], x[6]]
  return Likelihood(x_aux,"IH","no")
  
def Likelihood_c_Pee(x):
  # Tc0, tau_c, Rc, tau, Ta0=0, tau_a=0, Ma=0, Pee, toff_K, toff_I, toff_B
  x_aux=[x[0], x[1], x[2], x[3], 0, 0, 0, x[4], x[5], x[6], x[7]]
  return Likelihood(x_aux,"no","no")
 

############################################### Time-Integrated ##########################################

def param_profile(m,param,min_val,max_val,size=50):
    values=np.linspace(min_val,max_val,size)
    like,m_BFP=[],[]
    
    intial_values=np.array(m.values)
    m.fixed[param] = True

    for value in values:
        m.values=intial_values
        m.values[param]=value
        m.migrad()
        m_BFP.append(np.array(m.values))
        if not m.valid:
            print("Faield to converge in ",param," = ",value)
        like.append(m.fval)
    m.fixed[param] = False
    
    like=np.array(like)
    like=like-np.amin(like)
    
    return values,like,m_BFP
