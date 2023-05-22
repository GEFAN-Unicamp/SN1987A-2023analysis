import math
import numpy as np
from scipy.special import gamma
from flavor_conversion import * 

c = 1
h = 2*np.pi
D=50 #kpc
D = D * 100 * 3.086e19 # kpc to cm
Msun = 1.98848e30 * 1/(1.79e-30) # Sun Mass in MeV
Yn = 0.6	#Neutron fraction in the accreating material
mn = 939.565 #Neutron mass in MeV
erg_to_MeV=6.2415*10**5

######### Cooling #######################
def Tc(t, Tc0, tau_c):
  try:
    T = Tc0 * np.exp(-t/(4*tau_c))
  except RuntimeWarning:
    T = 10**(-20)
    print('Overflow in Tc function')
  return T
    
def ge_antnue(t, Enu, Tc0, tau_c):
	try:
	 	ge = Enu**2/(1 + np.exp(Enu/Tc(t, Tc0, tau_c)))
	except RuntimeWarning:
	  	print('Overflow in ge_antnue function, with Tc=',Tc(t, Tc0, tau_c))
	  	if type(Enu) == np.ndarray or type(Enu) == list:
	  		ge = np.zeros(len(Enu))
	  	else:
	  		ge = 0
	return ge
  
#cooling flux electron antineutrino
#Initially in MeV^2. Converting:
conv = (1e6 * 5.068e4)**2 * 1.519e15 # MeV^2 to 1/(cm² s eV)
def col_flux(t, Enu, Tc0,tau_c, Rc):
  #t[s], Enu[MeV], Tc0[MeV], tau_c[s], Rc[km]
  #if type(t) != np.ndarray:
  #  if t < 0:
      #return np.ones(len(Enu)) * 1e-20
  #    return 1e-20
  Rc = Rc * 1000 * 100 #km to cm
  flux = 1/(4*math.pi*D**2) * math.pi*c/(h*c)**3 * 4*math.pi * Rc**2 * ge_antnue(t, Enu, Tc0, tau_c) * conv # 1/(cm² s eV)
  return flux
  
 ######## Accretion ##############
def Ta(t, Ta0, Tc0, tau_a):
	Ti = Ta0
	Tf = 0.6*Tc0
	parm_m = 2
	Ta_t = Ti + np.sqrt((Tf - Ti)**2) * (t/tau_a)**parm_m
	return Ta_t
  
#time dependent function to guarantee continuity
#adm
def jk(t, tau_a):
  k_model = 2
  jk = np.exp(-(t/tau_a)**k_model)
  return jk

#neutron distribution function
#adm
def Nn(t, Ta0, Tc0, Ma, tau_a):
  Yn = 0.6
  mn = 939.565 #MeV
  Ma_Msun = Ma*Msun #Msun
  # Nn = Yn/mn * Ma_Msun * (Ta0/Ta(t, Ta0, Tc0, tau_a))**6 * jk(t, tau_a)/(1. + t/0.5) #with Ta(t)
  Nn = Yn/mn * Ma_Msun * jk(t, tau_a)/(1. + t/0.5) #with Ta cte
  return Nn

def Nn_2(t, Ta0, Tc0, Ma, tau_a):
  Yn = 0.6
  mn = 939.565 #MeV
  Ma_Msun = Ma*Msun #Msun
  Nn = Yn/mn * Ma_Msun * (Ta0/Ta(t, Ta0, Tc0, tau_a))**6 * jk(t, tau_a)/(1. + t/0.5) #with Ta(t)
  return Nn
 
 
#positron-neutron cross section
#cm^2
def cross_e_n(Enu):
  sigma_cm = 4.8e-44 * Enu**2 / (1 + Enu/260) #cm^2
  return sigma_cm


#positron energy
#MeV
def Epos_vec(Enu):
    Enu=np.array(Enu)
    Epos = np.where(Enu < 1.293, 0., (Enu - 1.293)/(1 - Enu/mn))
    return Epos

#positron spectrum
#MeV^2
def ge(t, Enu,Ta0, Tc0, tau_a):
  try:
    ge = Epos_vec(Enu)**2/(1 + np.exp(Epos_vec(Enu)/Ta0))
  except RuntimeWarning:
    ge = 0
  return ge

def ge_2(t, Enu,Ta0, Tc0, tau_a):
  try:
    ge = Epos_vec(Enu)**2/(1 + np.exp(Epos_vec(Enu)/Ta(t, Ta0, Tc0, tau_a))) #with Ta(t)
  except RuntimeWarning:
    ge = 0
  return ge
  
#accretion flux
#Initially in MeV^2. Converting:
conv = (1e6 * 5.068e4)**2 * 1.519e15 # MeV^2 to 1/(cm² s eV)
def acc_flux(t, Enu, Ta0, Tc0, tau_a, Ma):
#  if type(t) != np.ndarray:
#    if t < 0:
#      return np.ones(len(Enu)) * 1e-20
  flux = (8*math.pi*c)/(4*math.pi*D**2 * h**3 * c**3) *  Nn(t, Ta0, Tc0, Ma, tau_a) * cross_e_n(Enu) * ge(t, Enu, Ta0, Tc0, tau_a) * conv # 1/(cm² s eV)
  return flux

def acc_flux_2(t, Enu, Ta0, Tc0, tau_a, Ma):
  flux = (8*math.pi*c)/(4*math.pi*D**2 * h**3 * c**3) *  Nn_2(t, Ta0, Tc0, Ma, tau_a) * cross_e_n(Enu) * ge_2(t, Enu, Ta0, Tc0, tau_a) * conv # 1/(cm² s eV)
  return flux

#Garching fluence (time integrated flux)
def phi_Garching(E,E_0,L,alpha): 
    # E[MeV], E_0[MeV], L[10⁵²ergs], alpha
    L=L*(10**52)*erg_to_MeV
    N=((alpha+1)**(alpha+1))/(E_0*gamma(alpha+1))
    R=(L/E_0)*N*((E/E_0)**alpha)*np.exp((-1)*(alpha+1)*E/E_0)/(4*math.pi*D**2)
    return R #MeV⁻¹.cm⁻²


############################################################################################
def antinu_e_flux(t,Enu,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no"):
  #Antinu_e and Antinu_x initial fluxes
  #Antinu_x flux: only cooling with T_c0=tau*Tc0
  if acc=="yes":
    antinu_e_flux_0=acc_flux(t, Enu,Ta0, Tc0, tau_a, Ma) + (1 - jk(t, tau_a)) * col_flux(t - tau_a, Enu, Tc0, tau_c, Rc)
    antinu_x_flux_0=(1 - jk(t, tau_a)) * col_flux(t - tau_a, Enu, tau*Tc0, tau_c, Rc)  
  elif acc=="yes2":
    antinu_e_flux_0=acc_flux_2(t, Enu,Ta0, Tc0, tau_a, Ma) + (1 - jk(t, tau_a)) * col_flux(t - tau_a, Enu, Tc0, tau_c, Rc)
    antinu_x_flux_0=(1 - jk(t, tau_a)) * col_flux(t - tau_a, Enu, tau*Tc0, tau_c, Rc)
  elif acc=="yes_contemporaneous":
    antinu_e_flux_0=acc_flux(t, Enu,Ta0, Tc0, tau_a, Ma) + col_flux(t, Enu, Tc0, tau_c, Rc)
    antinu_x_flux_0=col_flux(t, Enu, tau*Tc0, tau_c, Rc)  
  elif acc=="no":
      antinu_e_flux_0=col_flux(t, Enu, Tc0, tau_c, Rc)
      antinu_x_flux_0=col_flux(t, Enu, tau*Tc0, tau_c, Rc)
  else:
    print('Invalid!')
    return 0

  Pee_aux=Pe_surv(Pee,mass_ord)

  antinu_e_flux= antinu_e_flux_0*Pee_aux+(1-Pee_aux)*antinu_x_flux_0
  return antinu_e_flux #MeV⁻¹.s⁻¹.cm⁻²


def antinu_e_flux_Garching(Enu,E_0_e,L_e,alpha_e,E_0_x,L_x,alpha_x,Pee,mass_ord,E_0_e_2=5,L_e_2=0,alpha_e_2=2.3,acc="no"):
  #Obs: same alpha for both antinu_e and antinu_x

  if acc=="yes":
    antinu_e_flux_0=phi_Garching(Enu,E_0_e,L_e,alpha_e)+phi_Garching(Enu,E_0_e_2,L_e_2,alpha_e_2)
    antinu_x_flux_0=phi_Garching(Enu,E_0_x,L_x,alpha_x)
  elif acc=="no":
    antinu_e_flux_0=phi_Garching(Enu,E_0_e,L_e,alpha_e)
    antinu_x_flux_0=phi_Garching(Enu,E_0_x,L_x,alpha_x)

  Pee_aux=Pe_surv(Pee,mass_ord)

  antinu_e_flux= antinu_e_flux_0*Pee_aux+(1-Pee_aux)*antinu_x_flux_0
  return antinu_e_flux #MeV⁻¹.cm⁻²
