import scipy.integrate as integrate
from flavor_conversion import*
import numpy as np
import math

Na=6.02*10**23 #Avogadro constant
H20_mol_mass=18.01528 #g/mol
C9H20_mol_mass=128.25#g/mol


#Dedin
#From prositron energy to neutrino energy
'''
def Enu_IB_vec(E_det,cos_theta): # E[MeV], cos(theta)
  E_det=np.array(E_det)
  cos_theta=np.array(cos_theta)
    
  delta=1.294 # Mass difference of a neutron and a proton [MeV]
  m_p=938.27208816 #Proton mass [MeV]
  m_pos=0.510998918 #Positron mass [MeV]
  E_pos_vec=E_det-m_pos #Subtracting the energy of the mass of the annihilated electron
    
  E_nu=[]
  #p_pos = np.where(E_pos_vec <=m_pos, 0.,np.sqrt((E_pos**2)-(m_pos**2)))
  #E_nu = np.where(p_pos <=0, 0.(E_pos+delta)/(1-((E_pos-(p_pos*cos_theta))/m_p)))
  for E_pos in E_pos_vec:
      if E_pos < m_pos:
          E_nu.append(0)
      else:
          p_pos=np.sqrt((E_pos**2)-(m_pos**2)) #positron momentum [MeV]
          E_nu_aux=(E_pos+delta)/(1-((E_pos-(p_pos*cos_theta))/m_p))
          if E_nu_aux< 1.806: #Energy Threshold [MeV] - Eq.(12.15) Giunti
            E_nu_aux=0
          E_nu.append(E_nu_aux)
  return E_nu
'''

def Enu_IB(E_det,cos_theta): # E[MeV], cos(theta)
  E_det=np.array(E_det)
  cos_theta=np.array(cos_theta)
   
  delta=1.294 # Mass difference of a neutron and a proton [MeV]
  m_p=938.27208816 #Proton mass [MeV]
  m_pos=0.510998918 #Positron mass [MeV]
  E_pos=E_det-m_pos #Subtracting the energy of the mass of the annihilated electron
  
  p_pos = np.where(E_pos <= m_pos, 0., np.sqrt((E_pos**2)-(m_pos**2))) #positron momentum [MeV]
  E_nu = np.where(p_pos <= 0, 0., (E_pos+delta)/(1-((E_pos-(p_pos*cos_theta))/m_p)))
  E_nu = np.where(E_nu <= 1.806, 0., E_nu) #Energy Threshold [MeV] - Eq.(12.15) Giunti
  
#   print(type(p_pos), type(E_nu))
#   if E_pos[0] < m_pos:
#     print(E_pos[0], m_pos, E_pos[0] < m_pos)
  
  return E_nu
#Enu_IB_vec=np.vectorize(Enu_IB)

'''
def Enu_IB(E_pos,cos_theta): # E[MeV], cos(theta)
  delta=1.294 #[MeV]
  m_p=938.27208816 #Proton mass [MeV]
  m_pos=0.510998918 #Positron mass [MeV]
  if E_pos < m_pos:
  	print(E_pos)
  	print("Warning: positron energy less than its mass!")
  	return 0
  p_pos=np.sqrt((E_pos**2)-(m_pos**2)) #positron momentum [MeV]
  E_nu=(E_pos+delta)/(1-((E_pos-(p_pos*cos_theta))/m_p))
  if E_nu< 1.806: #Energy Threshold [MeV] - Eq.(12.15) Giunti
    E_nu=0
  return E_nu
Enu_IB_vec=np.vectorize(Enu_IB)
'''

'''
#Marcos
#neutrino energy based on positron measured energy
def Eneutrino(Ee, cos_theta):
  mp = 938.272 #MeV
  me = 0.510999 #MeV
  ve = np.sqrt(1 - me**2/Ee**2)
  gamma = 1/np.sqrt(1 - ve**2)
  pe = gamma * me * ve
  try:
    Enu = (Ee + 1.294)/(1 - (Ee - pe * cos_theta)/mp)
  except IndexError:
    Enu = (Ee + 1.294)/(1 - (Ee - pe)/mp)
  return Enu
'''

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

'''
def dEnu_dEpos(Epos,cos_theta):
  delta_E=0.0001
  delta_f=Enu_IB_vec(Epos+delta_E/2,cos_theta)-Enu_IB_vec(Epos-delta_E/2,cos_theta)
  slope=delta_f/delta_E
  return slope
'''  
#angular IMB bias
def angular_bias(cos_theta):
  bias = 1 + 0.1*cos_theta
  return bias
  
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
  
#CC cross section
#cm^2
def cross_CC(Enu):
  sigmaCC = 9e-44 * Enu**2 #cm^2
  return sigmaCC

#ES cross section
#cm^2
def cross_ES(Enu):
  sigmaES = 9e-45 * Enu #cm^2
  return sigmaES
  
#The functions are from: https://www.sciencedirect.com/science/article/abs/pii/0370269388916516  
def eff(E,detec):
	if detec == 'K':
		eff=0.95*(1-np.exp(-1*((E-1.293)/9.3)**4))
	elif detec == 'I':
		eff=1-3*np.exp(-1*((E-1.293)/16))
		if type(eff) == np.ndarray:
			eff[eff<0]=0
		else:
			if eff<0:
				eff=0 	
	elif detec == 'B':
		eff=0.8*(1-np.exp(-1*((E-1.293)/10.17)**9))
	else:
		print("Not a detector option!")
		return 0
	return eff

###################################### Time-Dependent #################################################################

##### Kamiokande-II ##########
#COOLING + ACCRETION (d2N/dEdt) MeV-¹ s-¹
def dN_Kamiokande(t, E_det,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no"):
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
	dN_IB = np.where(Enu_detc_IB <=0,0.0, eff(E_det,'K')*dEnu_dEe(E_det,cos_theta)*Np*dsigma_dcos_vogel(Enu_detc_IB,cos_theta)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)*10**6) #eV-¹ s-¹ to MeV-¹ s-¹
  #dN_ES = Ne*cross_ES(E_det)*antinu_e_flux(t, E_det,Tc0,tau_c,Rc,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	return dN_IB


'''
def dN_Kamiokande(t, E_det,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no"):
	#http://arxiv.org/abs/astro-ph/0107260 - Lamb & Loredo		
	#Mass and number of targets
	#https://lss.fnal.gov/conf/C890928/p297.pdf
	fiducial_mass = 680e3 #kg
	fiducial_mass=2140 #ton
	fiducial_mass=fiducial_mass*10**6 #g
	Np_mol=2 #2 free protons (H) per molecule
	#Ne_mol=10 #electrons per molecule=number of protons (neutral material)
	Np=Np_mol*fiducial_mass*Na/H20_mol_mass #Number of free protons
	#Ne=Ne_mol*fiducial_mass*Na/H20_mol_mass #Number of free protons

	#Inverse beta decay event rate
	Enu_detc_IB=Enu_IB_vec(E_det,cos_theta)
	if Enu_detc_IB<=0:
		return 0
 	#OBS: 1/2 factor considering that the CS is isotropic over cos (theta)
	#dN_IB = eff(E_det,'K')*(1/2)*dEnu_dEe(E_det,cos_theta)*Np*cross_CC(Enu_detc_IB)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	dN_IB = eff(E_det,'K')*dEnu_dEe(E_det,cos_theta)*Np*dsigma_dcos_vogel(Enu_detc_IB,cos_theta)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)*10**6 #eV-¹ s-¹ to MeV-¹ s-¹
  #dN_ES = Ne*cross_ES(E_det)*antinu_e_flux(t, E_det,Tc0,tau_c,Rc,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	return dN_IB
  
dN_Kamiokande_vec=np.vectorize(dN_Kamiokande)
'''

##### IMB ##########
def dN_IMB(t, E_det,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no"):
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
	dN_IB = np.where(Enu_detc_IB <=0,0.0, eff(E_det,'I')*dEnu_dEe(E_det,cos_theta)*Np*dsigma_dcos_vogel(Enu_detc_IB,cos_theta)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)*10**6) #eV-¹ s-¹ to MeV-¹ s-¹
  #dN_ES = Ne*cross_ES(E_det)*antinu_e_flux(t, E_det,Tc0,tau_c,Rc,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	return dN_IB

'''
def dN_IMB(t, E_det,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no"):
  	#http://arxiv.org/abs/astro-ph/0107260 - Lamb & Loredo
  	#Mass and number of targets
  	#https://arxiv.org/pdf/1606.00665.pdf
	#https://www.researchgate.net/publication/23826760_A_study_of_atmospheric_neutrinos_with_the_IMB_detector
	fiducial_mass = 6800 #6800 tons
	fiducial_mass=fiducial_mass*10**6 #g
	Np_mol=2 #2 free protons (H) per molecule
	Np=Np_mol*fiducial_mass*Na/H20_mol_mass #Number of free protons
	
	#Inverse beta decay event rate
	Enu_detc_IB=Enu_IB_vec(E_det,cos_theta)
	if Enu_detc_IB<=0:
		return 0
	#OBS: 1/2 factor considering that the CS is isotropic over cos (theta)
	#dN_IB =  eff(E_det,'I')*angular_bias(cos_theta)*(1/2)*dEnu_dEe(E_det,cos_theta)*Np*cross_CC(Enu_detc_IB)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	dN_IB = eff(E_det,'I')*angular_bias(cos_theta)*dEnu_dEe(E_det,cos_theta)*Np*dsigma_dcos_vogel(Enu_detc_IB,cos_theta)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	return dN_IB
	
dN_IMB_vec=np.vectorize(dN_IMB)
'''
  	
##### Baksan ##########
def dN_Baksan(t, E_det,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no"):
	#http://arxiv.org/abs/astro-ph/0107260 - Lamb & Loredo
  	#https://www.sciencedirect.com/science/article/abs/pii/0370269388916516 - Baksan	
	#Mass and number of targets
	fiducial_mass=200 #ton
	fiducial_mass=fiducial_mass*10**6 #g
	Np_mol=20 #20 free protons (H) per C9H20 molecule
	Np=Np_mol*fiducial_mass*Na/C9H20_mol_mass #Number of free protons
	
	#Inverse beta decay event rate
	Enu_detc_IB=Enu_IB(E_det,cos_theta)
	dN_IB = np.where(Enu_detc_IB <=0,0.0, eff(E_det,'B')*dEnu_dEe(E_det,cos_theta)*Np*dsigma_dcos_vogel(Enu_detc_IB,cos_theta)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc)*10**6) #eV-¹ s-¹ to MeV-¹ s-¹
  #dN_ES = Ne*cross_ES(E_det)*antinu_e_flux(t, E_det,Tc0,tau_c,Rc,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	return dN_IB

'''
def dN_Baksan(t, E_det,cos_theta,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord="NH",acc="no"):
	#http://arxiv.org/abs/astro-ph/0107260 - Lamb & Loredo
  	#https://www.sciencedirect.com/science/article/abs/pii/0370269388916516 - Baksan	
	#Mass and number of targets
	fiducial_mass=200 #ton
	fiducial_mass=fiducial_mass*10**6 #g
	Np_mol=20 #20 free protons (H) per C9H20 molecule
	Np=Np_mol*fiducial_mass*Na/C9H20_mol_mass #Number of free protons
	
	#Inverse beta decay event rate
	Enu_detc_IB=Enu_IB_vec(E_det,cos_theta)
	if Enu_detc_IB<=0:
		return 0
 	#OBS: 1/2 factor considering that the CS is isotropic over cos (theta)
	#dN_IB = eff(E_det,'B')*(1/2)*dEnu_dEe(E_det,cos_theta)*Np*cross_CC(Enu_detc_IB)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	dN_IB = eff(E_det,'B')*dEnu_dEe(E_det,cos_theta)*Np*dsigma_dcos_vogel(Enu_detc_IB,cos_theta)* antinu_e_flux(t, Enu_detc_IB,Tc0,tau_c,Rc,tau,Ta0,tau_a,Ma,Pee,mass_ord,acc) * 10**6 #eV-¹ s-¹ to MeV-¹ s-¹
	return dN_IB
	
dN_Baksan_vec=np.vectorize(dN_Baksan)
'''
###################################### Time-Integrated #################################################################

#Kamiokande-II
def dN_Kamiokande_Garching(E_det,cos_theta,E_0_e,L_e,alpha_e,E_0_x,L_x,alpha_x,Pee,mass_ord):
  	#http://arxiv.org/abs/astro-ph/0107260 - Lamb & Loredo
  	#Mass and number of targets
  	#https://arxiv.org/pdf/1606.00665.pdf
	#https://www.researchgate.net/publication/23826760_A_study_of_atmospheric_neutrinos_with_the_IMB_detector
	#fiducial_mass = 3300 #3300 tons
	fiducial_mass=2.14*10**3 #ton
	fiducial_mass=fiducial_mass*10**6 #g
	Np_mol=2 #2 free protons (H) per molecule
	Np=Np_mol*fiducial_mass*Na/H20_mol_mass #Number of free protons

    #Inverse beta decay event rate
	Enu_detc_IB=Enu_IB(E_det,cos_theta)
	dN_IB = np.where(Enu_detc_IB <=0,0.0,eff(E_det,'K')*dEnu_dEe(E_det,cos_theta)*Np*cross_CC(Enu_detc_IB)* antinu_e_flux_Garching(Enu_detc_IB,E_0_e,L_e,alpha_e,E_0_x,L_x,alpha_x,Pee,mass_ord))
	return dN_IB #MeV⁻¹
#dN_Kamiokande_Garching_vec=np.vectorize(dN_Kamiokande_Garching)

##### IMB ##########
def dN_IMB_Garching(E_det,cos_theta,E_0_e,L_e,alpha_e,E_0_x,L_x,alpha_x,Pee,mass_ord):
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
	dN_IB = np.where(Enu_detc_IB <=0,0.0,eff(E_det,'I')*dEnu_dEe(E_det,cos_theta)*Np*cross_CC(Enu_detc_IB)* antinu_e_flux_Garching(Enu_detc_IB,E_0_e,L_e,alpha_e,E_0_x,L_x,alpha_x,Pee,mass_ord))
	return dN_IB #MeV⁻¹
#dN_IMB_Garching_vec=np.vectorize(dN_IMB_Garching)

##### Baksan ##########
def dN_Baksan_Garching(E_det,cos_theta,E_0_e,L_e,alpha_e,E_0_x,L_x,alpha_x,Pee,mass_ord):
	#http://arxiv.org/abs/astro-ph/0107260 - Lamb & Loredo
	#https://www.sciencedirect.com/science/article/abs/pii/0370269388916516 - Baksan	
	#Mass and number of targets
	fiducial_mass=200 #ton
	fiducial_mass=fiducial_mass*10**6 #g
	Np_mol=20 #20 free protons (H) per C9H20 molecule
	Np=Np_mol*fiducial_mass*Na/C9H20_mol_mass #Number of free protons
    
    #Inverse beta decay event rate
	Enu_detc_IB=Enu_IB(E_det,cos_theta)
	dN_IB = np.where(Enu_detc_IB <=0,0.0,eff(E_det,'B')*dEnu_dEe(E_det,cos_theta)*Np*cross_CC(Enu_detc_IB)* antinu_e_flux_Garching(Enu_detc_IB,E_0_e,L_e,alpha_e,E_0_x,L_x,alpha_x,Pee,mass_ord))
	return dN_IB #MeV⁻¹
#dN_Baksan_Garching_vec=np.vectorize(dN_Baksan_Garching)
