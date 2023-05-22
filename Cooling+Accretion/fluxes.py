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


  
#   #Flavor Conversion Mechanisms
#   if mass_ord=="NH":
#     Pee_aux= U[0][0]**2 #U^2_e1
#   elif mass_ord=="IH":
# #     Pf=P_f(Enu)
#     Pf=0
#     Pee_aux= (U[0][2]**2)*(1-Pf)+ (U[0][0]**2)*Pf #U^2_e3 * (1-Pf) + U^2_e1 * Pf
#   elif mass_ord=="no":
#     Pee_aux= Pee
#   elif mass_ord=="oqs":
#     #I have to set many functions to allow these parameters
#     gamma3_nat=1
#     gamma8_nat=1

#     D_kpc = 50
#     gamma3_kpc = gamma3_nat / (0.197e9 * 1e-15) * 3.086e19 # kpc-1
#     gamma8_kpc = gamma8_nat / (0.197e9 * 1e-15) * 3.086e19 # kpc-1
#     P11 = 1/3 + 1/2 * np.exp(-(gamma3_kpc + gamma8_kpc/3) * D_kpc) + 1/6 * np.exp(- gamma8_kpc * D_kpc)
#     #P22 = P11
#     #P33 = 1/3 + 2/3 * np.exp(-gamma8_kpc * D_kpc)
#     P21 = 1/3 - 1/2 * np.exp(-(gamma3_kpc + gamma8_kpc/3) * D_kpc) + 1/6 * np.exp(- gamma8_kpc * D_kpc)
#     P12 = P21
#     P31 = 1/3 - 1/3 * np.exp(-gamma8_kpc * D_kpc)
#     P13 = P31
#     #P32 = P31
#     #P23 = P31
#     #NH
#     #where Ueim ~ 1 because the initial matter effects
#     #Pee = P11 * U3[0][0]**2 + P12 * U3[0][1]**2 + P13 * U3[0][2]**2
#     Pee = U[0][0]**2 * P11 + U[0][1]**2 * P12 + U[0][2]**2 * P13
    
#   elif mass_ord=="E_swap":
#       Pee_aux_nubar=0
#       E_swap=Pee
#       Pee_aux=1-sig(Enu,E_swap,0.1)
        
#   elif mass_ord=="E_swap_inv":
#       Pee_aux_nubar=0
#       E_swap=Pee
#       Pee_aux=sig(Enu,E_swap,0.1)
        
#   else:
#     print('Invalid!')
#     return 0

#   antinu_e_flux= antinu_e_flux_0*Pee_aux+(1-Pee_aux)*antinu_x_flux_0
#   return antinu_e_flux #MeV⁻¹.s⁻¹.cm⁻²

# def antinu_e_flux_Garching(Enu,E_0_e,L_e,alpha_e,E_0_x,L_x,alpha_x,Pee,mass_ord,E_0_e_2=5,L_e_2=0,alpha_e_2=2.3,acc="no"):
#   #Obs: same alpha for both antinu_e and antinu_x

#   if acc=="yes":
#     antinu_e_flux_0=phi_Garching(Enu,E_0_e,L_e,alpha_e)+phi_Garching(Enu,E_0_e_2,L_e_2,alpha_e_2)
#     antinu_x_flux_0=phi_Garching(Enu,E_0_x,L_x,alpha_x)
#   elif acc=="no":
#     antinu_e_flux_0=phi_Garching(Enu,E_0_e,L_e,alpha_e)
#     antinu_x_flux_0=phi_Garching(Enu,E_0_x,L_x,alpha_x)

#   if mass_ord=="NH":
#     Pee_aux= U[0][0]**2 #U^2_e1
#   elif mass_ord=="IH":
# #     Pf=P_f(Enu)
#     Pf=0
#     Pee_aux= (U[0][2]**2)*(1-Pf)+ (U[0][0]**2)*Pf #U^2_e3 * (1-Pf) + U^2_e1 * Pf
#   elif mass_ord=="no":
#     Pee_aux= Pee
#   else:
#     print('Invalid!')
#     return 0

#   antinu_e_flux= antinu_e_flux_0*Pee_aux+(1-Pee_aux)*antinu_x_flux_0
#   return antinu_e_flux #MeV⁻¹.cm⁻²
