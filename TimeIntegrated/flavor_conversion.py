import numpy as np
#from fluxes import* 


# # Transition probability (not used)
# def P_f(E):
#   #if E==0:
#   #  return 0
#   Pf=np.exp(U[0][2]**2*((20/E)**(2/3))/(-3.5*10**(-5)))
#   return Pf
# #P_f_vec=np.vectorize(P_f)


def Pe_surv(Pee,mass_ord="NH"):

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

#Flavor Conversion Mechanisms
  if mass_ord=="NH":
    Pe_surv= U[0][0]**2 #U^2_e1
    return Pe_surv
  elif mass_ord=="IH":
    Pf=0
    Pe_surv= (U[0][2]**2)*(1-Pf)+ (U[0][0]**2)*Pf #U^2_e3 * (1-Pf) + U^2_e1 * Pf
    return Pe_surv
  elif mass_ord=="no":
    Pe_surv= Pee
    return Pe_surv
  else:
    print('Invalid!')
    return 0


