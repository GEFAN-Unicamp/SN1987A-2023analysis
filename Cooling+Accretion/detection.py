import scipy.integrate as integrate
#from flavor_conversion import*
from fluxes import * 
import numpy as np
import math
from scipy import special

Na=6.02*10**23 #Avogadro constant
H20_mol_mass=18.01528 #g/mol
C9H20_mol_mass=128.25#g/mol


#From prositron energy to neutrino energy
def Enu_IB(E_det,cos_theta): # E[MeV], cos(theta)
  E_det=np.array(E_det)
  cos_theta=np.array(cos_theta)
   
  delta=1.294 # Mass difference of a neutron and a proton [MeV]
  m_p=938.27208816 #Proton mass [MeV]
  m_pos=0.510998918 #Positron mass [MeV]
  E_pos=E_det-m_pos #Subtracting the energy of the mass of the annihilated electron
  
  p_pos = np.where(E_pos <= m_pos, 0., np.sqrt(abs((E_pos**2)-(m_pos**2)))) #positron momentum [MeV]
  E_nu = np.where(p_pos <= 0, 0., (E_pos+delta)/(1-((E_pos-(p_pos*cos_theta))/m_p)))
  E_nu = np.where(E_nu <= 1.806, 0., E_nu) #Energy Threshold [MeV] - Eq.(12.15) Giunti
  
#   print(type(p_pos), type(E_nu))
#   if E_pos[0] < m_pos:
#     print(E_pos[0], m_pos, E_pos[0] < m_pos)
  
  return E_nu

#derivative of neutrino energy dEnu/dEe
def dEnu_dEe(Ee, cos_theta):
  mp = 938.27208816 #Proton mass [MeV]
  me = 0.510998918 #Positron mass [MeV]
  pe = np.sqrt(Ee**2 - me**2)
  #ve = pe/Ee
  #ve = np.sqrt(1 - me**2/Ee**2)
  #gamma = 1/np.sqrt(1 - ve**2)
  #pe = gamma * me * ve
  #the derivative
  d = (1 + (1.294 + pe * cos_theta)/mp) / (1 - (Ee - pe*cos_theta)/mp)**2
  return d

# #angular IMB bias
# def angular_bias(cos_theta):
#   bias = 1 + 0.1*cos_theta
#   return bias
  
#differential cross section following (VOGEL and BEACON) (https://arxiv.org/pdf/hep-ph/9903554.pdf)
def dsigma_dcos_vogel(Enu, cos_theta):
  #Fermi constant
  Gf = 1.1663787e-11 #MeV-²
  #vector and axial-vector coupling constants
  f = 1
  g = 1.26
  mn = 939.56542052 #MeV  #PDG (2019)
  mp = 938.27208816 #MeV  #PDG (2019)
  me = 0.51099895000 #MeV  #PDG (2019)
  Delta = mn - mp
  Ee0 = Enu - Delta
  pe0 = np.sqrt(Ee0**2 - me**2)
  ve0 = pe0/Ee0
  M = (mp+mn)/2
  y_squared = (Delta**2 - me**2)/2
  Ee1 = Ee0 * (1 - Enu/M * (1 - ve0 * cos_theta)) - y_squared/M
  pe1 = np.sqrt(Ee1**2 - me**2)
  ve1 = pe1/Ee1
  cos_theta_C = 0.974
  Delta_R_inner = 0.024
  sigma0 = Gf**2 * cos_theta_C**2/math.pi * (1 + Delta_R_inner)
  f2 = 3.706
  Gamma = 2*(f + f2)*g * ((2*Ee0 + Delta) * (1 - ve0*cos_theta) - me**2/Ee0) + \
          (f**2 + g**2) * (Delta*(1 + ve0*cos_theta) + me**2/Ee0) + \
          (f**2 + 3*g**2) * ((Ee0 + Delta) * (1 - cos_theta/ve0) - Delta) + \
          (f**2 - g**2) * ((Ee0 + Delta) * (1 - cos_theta/ve0) - Delta) * ve0 * cos_theta
  
  hc_square = (0.197e-10)**2 #MeV^(-2) to cm^2
  dsigma_dcos = sigma0/2 * ((f**2 + 3*g**2) + (f**2 - g**2) * ve1 * cos_theta) * Ee1 * pe1 - sigma0/2 * Gamma/M * Ee0 * pe0
  dsigma_dcos = hc_square * dsigma_dcos
  return dsigma_dcos
  
# #CC cross section
# #cm^2
# def cross_CC(Enu):
#   sigmaCC = 9e-44 * Enu**2 #cm^2
#   return sigmaCC

#ES cross section
#cm^2
def cross_ES(Enu):
  sigmaES = 9e-45 * Enu #cm^2
  return sigmaES

######### Efficiency ##################
#Energy Uncertainty Function
def sigma_E(E,detec):
    if detec == 'K':
        err=1.27*np.sqrt(E/10)+1*(E/10)
    elif detec == 'I':
        err=3*np.sqrt(E/10)+0.4*(E/10)
    elif detec == 'B':
        err=0*np.sqrt(E/10)+2*(E/10)
    else:
        print("Not a detector option!")
        return 0
    return err #[MeV]

#Threshold Function
def g(E_e,E_min,detec):
    x= (E_e-E_min)/(np.sqrt(2)*sigma_E(E_e,detec))
    return (1+special.erf(x))/2
  
#Intrissic Efficiency
def eta(E_det,detec,model="Vissani"):
    if detec == 'K':
        if model=="our":        
            eff=0.95*(1-np.exp(-1*((E_det)/9.3)**4))
        elif model=="Vissani":
            eff=np.where(E_det<2.6,0,0.93*(1-(0.2/E_det)-(2.5/E_det)**2))
    elif detec == 'I':
        if model=="our":        
            eff=1-3*np.exp(-1*((E_det)/16))
            eff=np.where(eff<0,0,eff)
        elif model=="Vissani":
            c_IMB=[0.369,0,0,-6e-4,1e-4]
            E_IMB=15 #MeV
            eff=0
            for i in range(len(c_IMB)):
                eff=eff+c_IMB[i]*((E_det/E_IMB)-1)**(i+1)
            eff=np.where(E_det<E_IMB,0,eff)
            eff=np.where(eff>1,1,eff)
    elif detec == 'B':
        if model=="our":        
            eff=0.8*(1-np.exp(-1*((E_det)/10.17)**9))
        elif model=="Vissani":
            eff=np.where(E_det<0,0,1)
    else:
        print("Not a detector option!")
        return 0
    return eff

#Efficiency
def eff(E_det,detec,model="Vissani",typ="intrinsic",E_min=7.5):
#     E_nu=E+1.293
    if typ=="intrinsic":
        eff=eta(E_det,detec,model)
    elif typ=="total":
        if detec == 'K':
            eff=eta(E_det,detec,model)*g(E_det,E_min,detec)
        elif detec == 'I':
            eff=eta(E_det,detec,model)*g(E_det,E_min,detec)
        elif detec == 'B':
            eff=eta(E_det,detec,model)*g(E_det,E_min,detec)
        else:
            print("Not a detector option!")
            return 0
    else:
        print("Not a valid type of efficiency!")
        return 0
    return eff

###############################################################################################################
############################################### Rates in experiments ##########################################


##### Kamiokande-II ##########
#COOLING + ACCRETION (d2N/dEdt) MeV-¹ s-¹
def dN_Kamiokande(t, E_det,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no",eff_typ="intrinsic"):
	t_=np.array(t)
	E_det=np.array(E_det)
	cos_theta=np.array(cos_theta)
    
	#http://arxiv.org/abs/astro-ph/0107260 - Lamb & Loredo		
	#Mass and number of targets
	#https://lss.fnal.gov/conf/C890928/p297.pdf
	# fiducial_mass = 680e3 #kg
	fiducial_mass=2140 #ton
	fiducial_mass=fiducial_mass*10**6 #g
	Np_mol=2 #2 free protons (H) per molecule
	#Ne_mol=10 #electrons per molecule=number of protons (neutral material)
	Np=Np_mol*fiducial_mass*Na/H20_mol_mass #Number of free protons
	#Ne=Ne_mol*fiducial_mass*Na/H20_mol_mass #Number of free protons

	#Inverse beta decay event rate
	Enu_detc_IB=Enu_IB(E_det,cos_theta)
	dN_IB = np.where(Enu_detc_IB <=0,0.0, eff(E_det,'K',typ=eff_typ,E_min=4.5)*dEnu_dEe(E_det,cos_theta)*Np*dsigma_dcos_vogel(Enu_detc_IB,cos_theta)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)*10**6) #eV-¹ s-¹ to MeV-¹ s-¹
  #dN_ES = Ne*cross_ES(E_det)*antinu_e_flux(t, E_det,Tc0,tau_c,Rc,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	return dN_IB


##### IMB ##########
def dN_IMB(t, E_det,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no",eff_typ="intrinsic"):
  	#http://arxiv.org/abs/astro-ph/0107260 - Lamb & Loredo
  	#Mass and number of targets
  	#https://arxiv.org/pdf/1606.00665.pdf
	#https://www.researchgate.net/publication/23826760_A_study_of_atmospheric_neutrinos_with_the_IMB_detector
	fiducial_mass = 6800 #6800 tons
	fiducial_mass=fiducial_mass*10**6 #g
	Np_mol=2 #2 free protons (H) per molecule
	Np=Np_mol*fiducial_mass*Na/H20_mol_mass #Number of free protons
	
	#Inverse beta decay event rate
	Enu_detc_IB=Enu_IB(E_det,cos_theta)
	dN_IB = np.where(Enu_detc_IB <=0,0.0,eff(E_det,'I',typ=eff_typ,E_min=15)*dEnu_dEe(E_det,cos_theta)*Np*dsigma_dcos_vogel(Enu_detc_IB,cos_theta)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)*10**6) #eV-¹ s-¹ to MeV-¹ s-¹
  #dN_ES = Ne*cross_ES(E_det)*antinu_e_flux(t, E_det,Tc0,tau_c,Rc,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	return dN_IB

  	
##### Baksan ##########
def dN_Baksan(t, E_det,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no",eff_typ="intrinsic"):
	#http://arxiv.org/abs/astro-ph/0107260 - Lamb & Loredo
  	#https://www.sciencedirect.com/science/article/abs/pii/0370269388916516 - Baksan	
	#Mass and number of targets
	fiducial_mass=200 #ton
	fiducial_mass=fiducial_mass*10**6 #g
	Np_mol=20 #20 free protons (H) per C9H20 molecule
	Np=Np_mol*fiducial_mass*Na/C9H20_mol_mass #Number of free protons
	
	#Inverse beta decay event rate
	Enu_detc_IB=Enu_IB(E_det,cos_theta)
	dN_IB = np.where(Enu_detc_IB <=0,0.0, eff(E_det,'B',typ=eff_typ,E_min=10)*dEnu_dEe(E_det,cos_theta)*Np*dsigma_dcos_vogel(Enu_detc_IB,cos_theta)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)*10**6) #eV-¹ s-¹ to MeV-¹ s-¹
  #dN_ES = Ne*cross_ES(E_det)*antinu_e_flux(t, E_det,Tc0,tau_c,Rc,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	return dN_IB



###############################################################################################################

#Gaussian in energy and exp. uncertainty
def Gaussian_E(E, Edata, sigmaData):
  gauss = np.exp(-(Edata - E)**2./2./sigmaData**2)/(sigmaData*math.sqrt(2.*math.pi))
  return gauss
Gaussian_E_vec=np.vectorize(Gaussian_E)


def dN_Kamiokande_gaussian_vec(E, t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no") :
  gauss_R = dN_Kamiokande(t, E,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc) * Gaussian_E(E, Epos, delta_Epos)
  return gauss_R

def dN_IMB_gaussian_vec(E, t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no") :
  gauss_R = dN_IMB(t, E,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc) * Gaussian_E(E, Epos, delta_Epos)
  return gauss_R

def dN_Baksan_gaussian_vec(E, t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no") :
  gauss_R = dN_Baksan(t, E,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc) * Gaussian_E(E, Epos, delta_Epos)
  return gauss_R

def dN_Kamiokande_exp(t,Epos, delta_Epos, cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no") :
	E_max=Epos+5*delta_Epos
	E_min=Epos-5*delta_Epos
	if E_min<0.52:
		E_min=0.52
	E_vec=np.linspace(E_min,E_max,101)
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


