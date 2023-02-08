import numpy as np
from fluxes import* 

#PMNS
#Mixing angles
s12 = np.sqrt(0.307)
c12 = np.cos(np.arcsin(s12))
s23 = np.sqrt(0.545) #Normal order 
c23 = np.cos(np.arcsin(s23))
s13 = np.sqrt(2.18e-2)
c13 = np.cos(np.arcsin(s13))


#PMNS matrix
# Assumindo CP zero e neutrinos de Dirac
U23 = np.array([ [1,0,0],[0,c23,s23],[0,-s23,c23] ])
U13 = np.array([ [c13,0,s13],[0,1,0],[-s13,0,c13] ])
U12 = np.array([ [c12,s12,0],[-s12,c12,0],[0,0,1] ])
U = U23 @ U13 @ U12

#Flavor conversion
def P_f(E):
  #if E==0:
  #  return 0
  Pf=np.exp(U[0][2]**2*((20/E)**(2/3))/(-3.5*10**(-5)))
  return Pf
#P_f_vec=np.vectorize(P_f)

#Sigmoid function for the Spectral Split analysis
def sig(x,x_0=0,lamb=1):
    return 1/(1 + np.exp(-(x-x_0)/lamb))

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
  
  #Flavor Conversion Mechanisms
  if mass_ord=="NH":
    Pee_aux= U[0][0]**2 #U^2_e1
  elif mass_ord=="IH":
    Pf=P_f(Enu)
    Pee_aux= (U[0][2]**2)*(1-Pf)+ (U[0][0]**2)*Pf #U^2_e3 * (1-Pf) + U^2_e1 * Pf
  elif mass_ord=="no":
    Pee_aux= Pee
  elif mass_ord=="oqs":
    #I have to set many functions to allow these parameters
    gamma3_nat=1
    gamma8_nat=1

    D_kpc = 50
    gamma3_kpc = gamma3_nat / (0.197e9 * 1e-15) * 3.086e19 # kpc-1
    gamma8_kpc = gamma8_nat / (0.197e9 * 1e-15) * 3.086e19 # kpc-1
    P11 = 1/3 + 1/2 * np.exp(-(gamma3_kpc + gamma8_kpc/3) * D_kpc) + 1/6 * np.exp(- gamma8_kpc * D_kpc)
    #P22 = P11
    #P33 = 1/3 + 2/3 * np.exp(-gamma8_kpc * D_kpc)
    P21 = 1/3 - 1/2 * np.exp(-(gamma3_kpc + gamma8_kpc/3) * D_kpc) + 1/6 * np.exp(- gamma8_kpc * D_kpc)
    P12 = P21
    P31 = 1/3 - 1/3 * np.exp(-gamma8_kpc * D_kpc)
    P13 = P31
    #P32 = P31
    #P23 = P31
    #NH
    #where Ueim ~ 1 because the initial matter effects
    #Pee = P11 * U3[0][0]**2 + P12 * U3[0][1]**2 + P13 * U3[0][2]**2
    Pee = U[0][0]**2 * P11 + U[0][1]**2 * P12 + U[0][2]**2 * P13
    
  elif mass_ord=="E_swap":
      Pee_aux_nubar=0
      E_swap=Pee
      Pee_aux=1-sig(Enu,E_swap,0.1)
        
  elif mass_ord=="E_swap_inv":
      Pee_aux_nubar=0
      E_swap=Pee
      Pee_aux=sig(Enu,E_swap,0.1)
        
  else:
    print('Invalid!')
    return 0

  antinu_e_flux= antinu_e_flux_0*Pee_aux+(1-Pee_aux)*antinu_x_flux_0
  return antinu_e_flux #MeV⁻¹.s⁻¹.cm⁻²

def antinu_e_flux_Garching(Enu,E_0_e,L_e,alpha_e,E_0_x,L_x,alpha_x,Pee,mass_ord):
  #Obs: same alpha for both antinu_e and antinu_x
  antinu_e_flux_0=phi_Garching(Enu,E_0_e,L_e,alpha_e)
  antinu_x_flux_0=phi_Garching(Enu,E_0_x,L_x,alpha_x)

  if mass_ord=="NH":
    Pee_aux= U[0][0]**2 #U^2_e1
  elif mass_ord=="IH":
    Pf=P_f(Enu)
    Pee_aux= (U[0][2]**2)*(1-Pf)+ (U[0][0]**2)*Pf #U^2_e3 * (1-Pf) + U^2_e1 * Pf
  elif mass_ord=="no":
    Pee_aux= Pee
  else:
    print('Invalid!')
    return 0

  antinu_e_flux= antinu_e_flux_0*Pee_aux+(1-Pee_aux)*antinu_x_flux_0
  return antinu_e_flux #MeV⁻¹.cm⁻²
