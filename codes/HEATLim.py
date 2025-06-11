# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 11:24:20 2021

@author: Gisel Guzman-Echavarria

# =============================================================================
# Objetive model:
# =============================================================================
   
In this model, Survivability and Livability are assessed via the heat exchanges 
required and possible by the human body as explained in Vanos et al (2023), 
A physiological approach for assessing human survivability and liveability
to heat in a changing climate

Those heat exchanges are obtained using the heat balance equation based on a
 partitional calorimetry approach as stated in Cramer & Jay (2018).
 J Appl Physiol; CORP: Partitional Calorimetry v2.0

Disclaimer: The variables in this model follow the International System of Units (SI),
However there are some equations in which temperature is required in degrees Celsius.

The implementation of this model as it is assume steady-state conditions in environmental 
and physiological variables (Tskin and sweat rate).

Documentation of functions in this module follows https://www.softwaretestinghelp.com/python-docstring-tutorial/

 
# =============================================================================
# Input data in Cramer & Jay (2018):
# =============================================================================

Subject Characteristics - Corporal 
    - Mass   = Mass oh human body (kg)
    - Height = height of the human body (m)
    - AD     = Dubois-Dubois surface corporal area (m2)
    
    
    - Tsk_C  = Mean skin temperature (⁰C)	
    - Ar_AD  = Effective radiative area of the body (dimentionless)	
    
    - Emm_sk = Area weighted emissivity of the clothed body surface (dimentionless)		
    		
    - Icl	 =  Insulation value (Iclo)		
    - Re_cl	 = Evaporative heat transfer resistance of the clothing layer (m2·kPa·W-1)				
    
	- M - W, with W=0   Energy expenditure based on Look-up tables (W) as in
    Ainsworth, B. E., Haskell, W. L., Herrmann, S. D., Meckes, N., Bassett, D. R., Tudor-Locke, C., Greer, J. L., Vezina, J.,
    Whitt-Glover, M. C., & Leon, A. S. (2011). 2011 compendium of physical activities: A second update of codes and MET values.
    Medicine and Science in Sports and Exercise, 43(8), 1575–1581. https://doi.org/10.1249/MSS.0b013e31821ece12
    			
    - wmax_condition =  Maximum skin wettedness (dimentionless)
    - Max sweat rate (L/hr)
    

Environmental Characteristics - Ambiental 
    - Ta_C = Ambient temperature (⁰C)
    - Tr_C=	Mean radiant temperature (⁰C)
    - RH=		Relative humidity (%)
    - PB_kPa=  	Barometric pressure (kPa)
    - Av_ms=     Air velocity (m·s-1) 


Constants				
    σ	5.67E-08	#W·m-2·K-4	Stefan-Boltzmann constant.
    Lh_vap_water	= 2430	#J·g-1	Heat of vaporisation of water at 30⁰C. 
    LR	= 16.5	#K·kPa-1	Lewis Relation.
   
###############################################################################
# Disclaimer, not all the convertion factors listed here were used, it was preferred 
# to use MetPy package to move between the diferent metrics of humidity.
###############################################################################

"""
# =============================================================================
# load needed packages
# =============================================================================
import numpy as np
import pandas as pd

# =============================================================================
# Defining constants
# =============================================================================
σ = 5.67E-08	#W·m-2·K-4	Stefan-Boltzmann constant.
LR	= 16.5	#K·kPa-1	Lewis Relation.

Lh_vap =  2426	#J·g-1	Heat of vaporisation of Sweat heat at 30⁰C. 
#Gagnon, D., & Crandall, C. G. (2018). Sweating as a heat loss thermoeffector.
# Handbook of Clinical Neurology, 156, 211–232. https://doi.org/10.1016/B978-0-444-63912-7.00013-8
   
# =============================================================================
# Ambiental transformations, convertions and/or estimations
# =============================================================================
def T_Celsius_to_Kelvin(T_c:float) -> float:
    '''This function converts temperature from degrees Celsius to degrees Kelvin.
    
    Parameters
    ----------
    T_c : float
        Ambient temperature in degrees Celsius

    Returns
    ------
    T_K : float
        ambient temperature in degrees Kelvin'''
    
    T_K = T_c + 273.15
    return T_K

def Abs_Hum_from_TaC_PakPa(T_C:float,Pa_kPa:float) -> float:
    '''This function estimate absolute humidity from ambient temperature and 
    atmospheric pressure.
    
    Parameters
    ----------
    T_c  : float
        Ambient temperature in degrees Celsius
    Pa_kPa : float
        Atmospheric pressure in kPa
    
    Returns
    -------
    Abs_Hum: float
        Absolute Humidity in kg·m-3'''    
    Abs_Hum = 2.17 * (Pa_kPa/(T_C+273.15))
    return Abs_Hum

def Av_kmh_from_Av_ms(Av_ms:float) -> float:
    '''This function converts wind speed from kilometers per hour to meters per second.
    
    Parameters
    ----------
    Av_kmh : float
        Wind speed in kilometers per hour
    
    Returns
    -------
    Av_ms : float
        Wind speed in meters per second    
    '''    
    Av_kmh = Av_ms*3.6
    return Av_kmh


def Psa_kPa_from_TaC(T_C:float) -> float:
    '''This function estimates Saturated vapour pressure from ambient temperature.
    
    Parameters
    ----------
    T_c  : float
        Ambient temperature in degrees Celsius
    
    Returns
    -------
    Psa_kPa: float
        Saturated vapour pressure in kPa'''    
    
    Psa_kPa = np.exp(18.956-(4030.18/(T_C+235)))/10
    return Psa_kPa


def Pv_kPa_from_Psa_RH(Psa_kPa:float,RH:float) -> float:
    '''This function estimates vapour pressure in air mixed with water vapor,
        from vapour preesure and relative humidity.
    
    Parameters
    ----------
    Psa_kPa  : float
        Saturated vapour pressure in kPa
    RH  : float
        Relative humidity in percentage
    
    Returns
    -------
    Pv_kPa : float
        water vapour pressure (for air mixed with water vapor) in kPa'''    
    Pv_kPa = Psa_kPa*RH/100
    return Pv_kPa


def MET_to_MetabolicRate_W_Mass(METS:float, Mass:float) -> float:
    '''The conversion factor from watts to metabolic equivalents (METs) and 
    vice versa is taken as from MET definition as the amount of oxygen consumed
    while sitting at rest as 3.5 ml O2/kg.min (Jetté et al., 1990), 
    while assuming a fixed energy yield of 21 kJ/l O2, 
    (i.e. respiratory exchange ratio of ~0.97). 
    The combination of both factors leads a conversion factor of 1 MET=1.225 W/kg as a result.
    
    Parameters
    ----------
    M_MET : float
        Energy expenditure in METS
    Mass : float
         Body mass in kg 
       
    Returns
    -------
    M_W : float
        Energy expenditure in Watts
    ''' 
    Watts = METS*1.225*Mass
    return Watts


def MetabolicRate_W_to_MET_Mass(Watts:float, Mass:float) -> float:
    '''The conversion factor from watts to metabolic equivalents (METs) and 
    vice versa is taken as from MET definition as the amount of oxygen consumed
    while sitting at rest as 3.5 ml O2/kg.min (Jetté et al., 1990), 
    while assuming a fixed energy yield of 21 kJ/l O2, 
    (i.e. respiratory exchange ratio of ~0.97). 
    The combination of both factors leads a conversion factor of 1 MET=1.225 W/kg as a result.
        
    Parameters
    ----------
    M_WAT : float
       Energy expenditure in METS
    Mass : float
       Body mass in kg 
        
    Returns
    -------
    M_MET : float
        Energy expenditure in Watts'''    

    METS = Watts/(1.225*Mass)
    return METS

# =============================================================================
# reading input data
# =============================================================================
def read_personal_profiles(path_profile:str) -> dict:
    '''This functions reads the personal profile of the person to run the model.
    
    Parameters
    ----------
    path_profile : string
        Path+name of the text file in which the personal profiles is located in the computer
    
    Returns
    -------
    profile : dict
        Information about the personal profile to run the model'''    
    dtype_profile={'id':str, 'name':'str', 'name_file':'str', 
                   'Mass':float, 'Height':'float32', 'AD':'float32', 
                   'Tsk_C':'float32', 'A_eff':'float32', 'Emm_sk':'float32', 
                   'Icl':'float32', 'Re_cl':'float32','M':'float32', 'W':'float32', 
                   'wmax_condition':'str', 'wmax_rate':'float32'} 
    profile = pd.read_csv(path_profile,header = None,index_col = 0,dtype=dtype_profile, sep = '\t',na_values= -9999,names = ['values'])
    profile = profile.to_dict()['values']
    return profile

# =============================================================================
# personal transformations, estimations and/or corrections 
# =============================================================================


def AD_from_mass_height(Mass:float,Height:float) -> float:  
    '''This function returns Dubois-Dubois surface corporal area.
    
    Parameters
    ----------
    Mass : float
        Mass in kg
    Height: float
        Height in meters
    
    Returns
    -------
    ad : float
        Corporal surface area based on Dubois-Dubois equation'''    
    ad =0.202*((Mass**0.425)*((Height**0.725)))
    return ad

    
# =============================================================================
# Dry Heat Loss from the Skin (C+R)		
# =============================================================================

def fcl_from_Icl(Icl:float) -> float:
    '''This function estimated the clothing area factor (dimensionless).
    
    Parameters
    ----------
    Icl : float 
        Insulation value (Iclo)	
    
    Returns
    -------
    fcl : float
        Clothing area factor'''    
    fcl = 1 + (0.31*Icl)
    return fcl

def hc_cof_from_Av(Av_ms:float) -> float:
    '''This function estimate convective heat transfer based on:
    
    Parsons, K. (2014). Human Thermal Environments (Third edit). CRC Press. 
    https://doi.org/10.1201/b16750 (page 51)
    
    Parameters
    ----------
    Av_ms : float
        Wind speed in meters per second
    Returns
    -------
    hc_cof : float
        Convective heat transfer coefficient in W/(m2.C)'''    
    hc_cof = np.where(Av_ms < 0.2, 3.61,((Av_ms**0.6)*8.3))
    return hc_cof

def Rcl_from_Icl(Icl:float) -> float:
    '''This function estimates dry heat transfer resistance of clothing
     from “clo” units (1 clo = 0.155 m2.°C/W).
    
    Parameters
    ----------
    Icl : float
        Insulation value (Iclo)	
    Returns
    -------
    Rcl : float
        Dry heat transfer resistance of clothing in m2.°C/W'''    
    #Intrinsic clothing insulation
    Rcl = 0.155*Icl
    return Rcl
        
def hr_cof_from_radiant_features(Tr_C:float,Tsk_C:float,Emm_sk:float,A_eff:float) -> float:
    '''This function estimates radiative heat transfer coefficient.
    
    Parameters
    ----------
    Tr_C : float
        Mean radiant temperature in degress Celsius
    Tsk_C : float
        Mean skin temperature in degreess Celsius
    Emm_sk : float
        rea weighted emissivity of the clothed body surface (dimentionless)	
    A_eff : float
        Effective radiative area of the body (dimentionless).

    Returns
    -------
    hr_coef : float
        Radiative heat transfer coefficient W/m2·K'''    
    Boltzmann = 5.67*10**-8 #W·m-2·K-4	Stefan-Boltzmann constant.
    hr_cof = 4*Emm_sk*Boltzmann*A_eff*((273.2+(Tsk_C+Tr_C)/2)**3)

    return hr_cof
    

def h_coef_from_hc_hr(hc_cof:float,hr_cof:float) -> float:
    '''This function estimates combined convective heat transfer coefficient.
 
    Parameters
    ----------
    hc_coef : float
        Convective heat transfer coefficient in W·m-2·K-1
    hr_coef : float
        Radiative heat transfer coefficient in W·m-2·K-1   
        
    Returns
    h_cof : float
        combined Heat transfer coefficient in W·m-2·K-1
    -------
    : '''    
    h_cof = hc_cof + hr_cof
    return h_cof

def to_from_hr_tr_hc_ta(hr_cof:float,mrt_C:float,hc_cof:float,Ta_C:float) -> float:
    '''This function estimates operative temperature from the radiative heat transfer coefficient,
    mean radiant temperature, convective heat transfer coefficient and air temperature.
    
    
    Parameters
    ----------
    hr_coef : float
        Radiative heat transfer coefficient in W·m-2·K-1
    mrt_C: float
       Mean radiant temperature in degrees Celsius
    hr_coef : float
        Radiative heat transfer coefficient in W·m-2·K-1
    Ta_C: float
        Ambient temperature in degrees Celsius
        
    Returns
    -------
    to_C : float
        Operative temperature in degrees Celsius
        '''    
    to_C=((hr_cof*mrt_C)+(hc_cof*Ta_C))/(hr_cof+hc_cof)
    return to_C

def Dry_Heat_Loss_c_plus_r(Tsk_C:float,to_C:float,Icl:float,h_cof:float,Ad:float) -> float:
    '''This function  estimates combined dry heat loss via convection and radiation.
    
    Parameters
    ----------
    Tsk_C : float
        Mean skin temperature in degrees Celsius
    to_C : float
        Operative temperature in degrees Celsius
    Icl:float
        Insulation clothing value  in CLO	
    h_cof:float
        combined convective heat transfer coefficient in W·m-2·K-1
    Ad : float
        Corporal surface area in m2
        
    Returns
    -------
    Dry_Heat_Loss: float
        combined dry heat loss via convection and radiation in W'''    
    fcl	= fcl_from_Icl(Icl)	#ND	Clothing area factor. Ratio of clothed body surface to nude body surface.
    Rcl	= Rcl_from_Icl(Icl)	#m2·⁰C·W-1	Intrinsic clothing insulation
    
    
    Dry_Heat_Loss_AD = (Tsk_C-to_C)/(Rcl+(1/(h_cof*fcl)))
    Dry_Heat_Loss = Dry_Heat_Loss_AD*Ad #Convertion from w.M-2 to W
    return Dry_Heat_Loss


# =============================================================================
# Heat Loss from Respiration (Cres+Eres)	
# =============================================================================

	
def Cres_from_M_Ta(M:float,Ta_C:float,Ad:float) -> float:
    '''This function estimates the respiratory heat loss via convection, using
    ASHRAE 1997 in W/m2 units, then we need to multiply by surface area to
    obtain the heat flux in W.  
    
    Parameters
    ----------
    M : float
        Rate of metabolic energy expenditure in W
    Ta_C: float
        Ambient temperature in degrees Celsius
    Ad : float
        Corporal surface area in m2
        
    Returns
    -------
    Cres : float 
        Respiratory heat loss via convection in W'''    
    Cres = 0.0014*M*(34-Ta_C)*Ad
    return Cres

def Eres_from_M_Pa(M:float,Pa_kPa:float,Ad:float) -> float:
    '''This function estimates Latent respiratory heat loss.
    ASHRAE 1997 in W/m2 units, then multiply by surface area to obtain
    the heat flux in W.
    
    Parameters
    ----------
    M : float
        Rate of metabolic energy expenditure in W
    Pv_kPa : float
        water vapour pressure (for air mixed with water vapor) in kPa
    Ad : float
        Corporal surface area in m2
    Returns
    -------
    Eres : float
        Respiratory heat loss via evaporation in W'''    
    Eres = 0.0173*M*(5.86618428-Pa_kPa)*Ad
    
    return Eres


# =============================================================================
# Estimation of Evaporation required and ambiental possible
# =============================================================================


def Ereq_from_HeatFluxes(M:float,W:float,Dry_Heat_Loss:float,Cres:float,Eres:float) -> float:
   
    '''This function estimates the amount of evaporative heat loss required for heat balance.
    
    Ereq  = (M - Wk) - (C+R) - (Cres + Eres)	
    Ereq = Hprod - Dry_Heat_Loss - CEplus_res
    
    Parameters
    ----------
    M : float
        If W is correspond to internal heat production in W, otherwise to metabolic rate in W
    W : float
        External work done for human body in W
    Dry_Heat_Loss: float
        Heat dry heat transfer by radiation and convection trough the skin in W
    Cres : float
        Dry respiratory heat loss by convection in W
    Eres : float
        Latent respiratory heat loss in W
    
    Returns
    -------
    Ereq : float
        Rate of evaporation  in watts required for heat balance to whole body'''    
    
    Hprod = M - W
    CEplus_res = Cres + Eres
    Ereq = Hprod  - Dry_Heat_Loss - CEplus_res
    
        
    return Ereq

def he_cof(hc_cof:float) -> float:
    '''This function estimates the evaporative heat transfer coefficient using the
    Lewis relation (16.5 K/kPa).
    
    Parameters
    ----------
    hc_cof : float
        Convective heat transfer coefficient in W/m2·K
    Returns
    -------
    he_cof : float
        Evaporative heat transfer coefficient in W/m2·kPa'''    
    #As in partitional calorimetry model excel spreadsheet from Ollie Jay
    
    he_cof =hc_cof*LR
    return he_cof

def wmax(person_condition:str) -> float: 
    '''This function provide the maximum skin wettedness depending on the chareacteristic
    set in the personal profile file:
    
    ISO:
    Unacclimated = 0.85
    Acclimatied = 1.00
    
    Ravanelli et al. MSSE (2018):
    Untrained & Unacclimated = 0.72
    Trained & Unacclimated = 0.84
    Trained & Acclimated = 0.95
    
    Morris 2021
    0.85 for the YNG model (Candas et al., 1979a);
    0.65 for the OLD model
    
    NOTE: This is a factor to be improved in the future once there is more data
    arounf from thermal physiologist.
        
    Parameters
    ----------
    person condition : str
        Describe if the person is acclimatized or not. Also if have sweating impairments or not.
    
    Returns
    -------
    wmax : float
        Maximum skin wettednes'''    
    
    if   person_condition == 'Unacclimated':  wmax = 0.85
    elif person_condition == 'fully acclimated':   wmax = 1
    elif person_condition == 'Untrained & Unacclimated':  wmax = 0.72
    elif person_condition == 'Trained & Unacclimated':  wmax = 0.84
    elif person_condition == 'Trained & Acclimated':  wmax = 0.95
    elif person_condition == 'YNG_Morris_2021':  wmax = 0.85
    elif person_condition == 'OLD_Morris_2021':  wmax = 0.65
    else: 
        print('Invalid "person_condition"')
        
    return wmax

def Emax_env(Psk_s:float,Pv_kPa:float,Re_cl:float,he_cof:float,Icl:float,AD:float) -> float:
    '''This function estimates the maximum evaporative heat loss for a given thermal environment 
    and clothing, also known as the biophysical evaporative heat loss.
    
    Parameters
    ----------
    Psk_s : float
        Vapor pressure at the skin surface while saturated with sweat in kPa
    Pv_kPa : float
        Ambient vapour pressure in kPa
    Re_cl : float
        Evaporative resistance of clothing in m2·kPa·W-1
    he_cof : float
        Evaporative heat transfer coefficient in W·m-2·kPa-1
    Icl : float
        
    AD : float
        Dubois surface corporal area in m2
    
    Returns
    -------
    Emax_env: float
        Biophysical evaporative heat loss (cause by ambient environment and the
        clothes people wear) in Watts
    '''
    
    fcl	= fcl_from_Icl(Icl)	#Ratio of clothed body surface to nude body surface (dimensionless)

    Emax_env_AD = (Psk_s-Pv_kPa)/(Re_cl+(1/(he_cof*fcl)))            
    Emax_env = Emax_env_AD*AD 
    return Emax_env

def Emax_wettedness(wmax:float,Psk_s:float,Pv_kPa:float,Re_cl:float,he_cof:float,Icl:float,AD:float) -> float:
    '''This function estimates the biophysical evaporative heat loss accounting as additional 
    constraing the capacity to physiologically wet the skin and thus distribute sweat across 
    the skin surface in humid environments. This restriction is applied using the maximum skin wettedness.
    
    Parameters
    ----------
    wmax : float
        Maximum or critical skin wettedness (dimensionless)
    Psk_s : float
        Vapor pressure at the skin surface while saturated with sweat in kPa
    Pv_kPa : float
        Ambient vapour pressure in kPa
    Re_cl : float
        Evaporative resistance of clothing in m2·kPa·W-1
    he_cof : float
        Evaporative heat transfer coefficient in W·m-2·kPa-1
    Icl : float
        
    AD : float
        Dubois surface corporal area in m2
    
    Returns
    -------
    Emax_wettedness: float
        Biophysical evaporative heat loss modified by the constrain of max wettedness in Watts
    '''
    
    fcl	= fcl_from_Icl(Icl)	#Ratio of clothed body surface to nude body surface (dimensionless)

    Emax_wettedness_AD = wmax*(Psk_s-Pv_kPa)/(Re_cl+(1/(he_cof*fcl)))            
    Emax_wettedness = Emax_wettedness_AD*AD 
    return Emax_wettedness


    
def Emax_sweat_rate(Smax:float,Lh_vap:float, density:float, r:float) -> float:
    '''This function estimates the evaporative heat loss after the evaporation 
    of the maximum volume of sweat that people can release based on the 
    maximum hourly sweat rate (Smax) and accounting sweat efficiency. (r)
    
    In this equation 3.6 is the conversion factor to account the time from hours 
    to seconds and the volume of sweat to mass, assuming that sweat density is similar like water
    
    This metric implies a constant sweat rate and a completely water replenish.
    
    Parameters
    ----------
    Smax : float
        Maximum sweat rate for a given personal profile in L/h
    Lh_vap: float
        Heat latent of vaporization of sweat in J·g-1 or the amount of energy in form of 
        enthalpy that is add to the air when sweat evaporate
    density: float
        Density of sweat, assumed here as 
    r : float
        Sweating efficiency (dimentionless)
    
    Returns
    -------
    Emax_sweat_rate : float
        Maximun evaporative heat loss in W linked with the sweat evaporation of the maximum sweat rate'''    
    Emax = ((Smax*Lh_vap*density)/3.6)*r
    return Emax
# =============================================================================
# Body capacity to sweat
# =============================================================================

#the body is capable? how stress the body is?

def wreq_HSI_skin_wettedness(Ereq:float,Emax_env:float) -> float:
    '''This function estimates biophysical skin wettedness required (ω_req) for heat balance.
    
    skin wettedness: the proportion of the skin surface saturated with sweat. 
    Notice that there is a maximum value of skin wettedness that determine when heat is compensable.
    
    Parameters
    ----------
    Ereq : float
        Evaporative heat rate required for heat balance  W
    Emax_env : float
        Biophysical evaporative heat loss (environment + clothing) in W
    
    Returns
    -------
    wreq: float
        Biophysical skin wettedness required for heat balance'''    

    wreq = Ereq/Emax_env	
    return wreq		
	

def Sweating_efficiency_r(wreq:float) -> float:
    '''This function estimates the sweating efficiency based on the skin wettedness.
    That value can be taken as the proportion of sweat produced that is not dripped off
    the body and evaporated from the skin surface, thus contributing to evaporative heat loss.
    A minimum of sweating efficiency was set at 0.55 based on (Candas et al., 1979a).
    
    Note that the minimum value is set to 0.5 at wreq = 1, for future applications in cold conditions
    ask if there is a need of truncate the r to 1.
    
    Parameters
    ----------
    wreq: float
        Skin wettedness also know as Heat Stress Index
    
    Returns
    -------
    r : float 
        Sweating efficiency (dimentionless)'''    

    r = np.where(wreq<1, 1-((wreq**2)/2),0.5)
    r = np.where(r>1,1,r)
    return r

def Sreq(Ereq:float,r:float,Lh_vap:float) -> float:
    '''This functions estimates the required sweat rate to maintain heat balance
    and therefore compensate an imposed thermal load.
    
    In this equation 3.6 is the conversion factor to account the time from hours 
    to seconds and the volume of sweat to mass, assuming that sweat density like
    water density.

    Parameters
    ----------
    Ereq : float
        Rate of evaporation required for heat balance to whole body
    Lh_vap: float
        heat latent of vaporization of sweat in J·g-1 or the amount of energy in form of 
        enthalpy that is add to the air when sweat evaporate
    r : float
        Sweating efficiency (dimentionless)
    
    Returns
    -------
    Sreq : float 
        Required sweat rate L·h-1'''    
    Sreq =((Ereq/r)/Lh_vap)*3.6
    return Sreq

# =============================================================================
# Survivability and livability 
# =============================================================================

def Ssurvive_from_Exposure_time(Exp_time):
    ''' This function estimates the critical rate of heat storage (Ssurv) before inevitable
    heat stroke death during rest for exposure times of 6 and 3 hours.
    Here is assumed to have a linear increase of 6°C (Starting temperature at 37 °C)
    and a Human body-specific heat capacity of C_p= 2.98 kJ/(Kg.°C). 
    
    Cp: Xiaojiang Xu, Timothy P. Rioux & Michael P. Castellani (2022) The specific heat of the human 
    body is lower than previously believed: The Journal Temperature toolbox, Temperature, 
    DOI: 10.1080/23328940.2022.2088034

    
    '''
    valid = {1, 3, 6}
    if Exp_time not in valid:
        raise ValueError("results: status must be one of %r." % valid)
    #Estimation of constant power along exposure time
    DQ = 2.98*6  #CP.dT (change of temperature of 6 degrees)
    Ssurvive = (DQ * 1000)/(Exp_time*3600) #1000 is conversion from KJ to J, and 3600 conversion seg to hours
    
    return Ssurvive

def Survivability(Exp_time:float,Ereq:float,Emax_wettedness:float,Emax_sweat:float,Sreq:float,Smax:float,Mass:float) -> tuple:
    '''This function asses Survibability on humans before inevitable
    heat stroke death during rest for exposure times of 6 and 3 hours
    in a given thermal environment. In the fuction Ssurvive_from_Exposure_time of this 
    very same module is the estimation of the critical rate of heat storage 
    (Ssurv) before death.
    
    ** This assesment assume People will survive, even if they cannot thermally
    compensate the environment and sometimes the heat storage can be greater than zero,
    while the body core temperature does not surpass 43°C
        
    ** See all the environmental and physiological assumptions in the main paper.
    
    After ends the decision-making process of the algorithm also is assigned a flag
    that categorized the survivability type according to the physiological 
    constraints imposed in this model (See in Supplemental material, Figure S1 Model workflows for (a)
        survivability and (b) liveability) 
        
    Thus, the algorithm determines survivability (as a dichotomous variable: yes/no) 
    and assigns outcomes based on combined environmental and physiological restrictions.
    
    Based on this framework, a person will (Notice the numbers represent the survivability 
                                            zones in the output flag_survivability): 
    
    1.	survive while remaining within sweating limits
    2.	survive despite exceeding sweating limits
    3.	not survive because the environment restricts heat loss too much (in high humidity)
    4.	not survive because the required sweat rate is not possible (in low humidity)
    5.	not survive due to both critical environmental heat loss restrictions (3rd argument) and not possible sweat rate to dissipate heat (4th argument).

    Parameters
    ----------
    Exp_time : float
        Exposure time in hours. Valid values are: 1H, 3H or 6H
    Ereq : float
        Rate of evaporation in Watts required for heat balance to whole body    
    Emax_wettedness : float
        Biophysical evaporative heat loss modified by the constrain of max wettedness in Watts
    Emax_sweat : float
        Maximun evaporative heat loss in Watts linked with the sweat evaporation of the maximum sweat rate    
    Sreq : float
        Required sweat rate in L·h-1        
    Smax : float
        Maximum sweat rate for a personal profile in L·h-1
    Mass:
        Mass of a person in kg
    
    Returns
    -------
    survivability : boolean
        Could a person with a given personal profile exposed to a certain survive a given thermal environment? yes, no
    flag_survivability : int
        Flag that indicates the survivability type according to the physiological constraints imposed in this model
    
    '''    
    valid = {1, 3, 6}
    if Exp_time not in valid:
        raise ValueError("results: status must be one of %r." % valid)
        
    Ssurvive = Ssurvive_from_Exposure_time(Exp_time)

    #For the explanation of the workflow or criteria followed here, also See
    #Figure S1 (Model workflows for (a) survivability and (b) liveability) in the suplemental material.      
    #Criteria 1:    
    condition1 = np.where((Ereq - Emax_wettedness) <=  (Ssurvive*Mass), True,False) #(W/kg) multiply by Mass to convert in W

    #Criteria 2.1:
    condition21 = np.where(Sreq <=  Smax, True,False) #Is sweat enough?
        
    #Criteria 2.2:
    condition22 = np.where((Ereq - Emax_sweat) <=  (Ssurvive*Mass), True,False)#(W/kg) multiply by Mass to convert in W

    #Assesing survivability:
    survivability = np.zeros(condition1.shape, dtype = bool) #by default people do not survive
    flag_survivability = np.zeros(condition1.shape)*np.nan
     
    
    #survivability[condition1 == False] = False
    flag_survivability[condition1 == False] = 3
    
    aux1 = np.logical_and(condition1, condition21) #if condition 1 and 2.1 are true
    survivability[aux1 == True] = True
    flag_survivability[aux1 == True] = 1
 
    aux2 = np.logical_and(condition1, ~condition21) #if condition is true 1 and 2.1 is false
    aux3 = np.logical_and(aux2, condition22) #Emax_sweat from Smax is enough
    survivability[aux3 == True] = True
    flag_survivability[aux3 == True] =2

    aux4 = np.logical_and(aux2, ~condition22) #Emax_sweat from Smax is not enough
    #survivability[aux4 == False] = False
    flag_survivability[aux4 == True] = 4
    
    aux5 = np.logical_and(~condition1, ~condition22)
    flag_survivability[aux5 == True] = 5

    return survivability, flag_survivability



def livability_Mmax(survivability: bool,Ereq:float,Emax_wettedness:float,Emax_sweat:float,M_rest:float, Mass:float) -> tuple:
    '''Liveability is the maximum metabolic rate (Mmax) that can be generated before S≥0, 
    or sustained compensable heat stress, whith M = Hprod. The Mmax value indicates the 
    sustained activity levels (intensity but not duration) possible without unchecked 
    rises in Tcore (i.e., uncompensable heat stress) within a given steady-state 
    thermal environment. 
    
    After checking that heat stress is compensable (Ereq≤Emaxlim), then Mmax is estimated 
    as follows.
    
        M_max  = Emax_lim  -  H_loss   
        
    Notice that:
        - Emax_lim is the minimum value among Emax_wettedness and Emax_sweat
        - H_loss can be represented either by (Hdry +Cres +Eres) or (Ereq + Hprod). 
    
    Also if Mmax is cero or less than cero means the thermal load is non-compensable, then that person 
    can survive but is not able to live (no activity is possible without storage heat internally).
    However, there is a limit and that is why the survivability assesment is accounted for this variable.

    
    Parameters
    ----------
    survivability : bool
        Boolean variable indicating if a person can survive in a given thermal environment
    Ereq : float
        Rate of evaporation in Watts required for heat balance to whole body    
    Emax_wettedness : float
        Biophysical evaporative heat loss modified by the constrain of max wettedness in Watts
    Emax_sweat : float
        Maximun evaporative heat loss in Watts linked with the sweat evaporation of the maximum sweat rate    
    M_rest : float
        Metabolic energetic expenditure while people is resting in W (here M_rest = than Hprod)
    Mass : float
        Body mass in kg 
    Returns
    -------
    Mmax : float  
        Maximum metabolic rate in W that a person can tolerate after surviving, 
        that reflects the range of activities could be performed
        
    mask_non_livable:
        Flag to indicate if this environmental conditions leads to a survivable but not
        liveable condition (no activity is possible without storage heat internally).
        '''   
    
    Emax_constrain = np.min([Emax_wettedness,Emax_sweat],axis=0)
    compensability = np.where(Ereq < Emax_constrain, True,False) 
    
    Mmax = Emax_constrain - Ereq + M_rest
    
    Mmax[~compensability] = np.nan
    non_livable = np.logical_and(survivability, ~compensability)
    Mmax_MET = MetabolicRate_W_to_MET_Mass(Mmax,Mass) #Convert to METs
    Mmax_capped = np.where(Mmax_MET > 10, 10, Mmax_MET) #if Mmax is greather than 10 METs
    return Mmax_capped, non_livable


# =============================================================================
# Lines extraction limits 
# =============================================================================

def SurvivabilityLines_from_SurvivabilityMatrix(Survivability:bool,humidity:float,type_humidity:str,temperature:float,tw:float) -> float:
    ''' This function extract the survivability limits from the survibability 
    results of the model (yes/no matrix), given the assesment results 
    from a temperature and humidity range.
    
    This works from assesments based on arrays of the temperature and humidity matrix, with
    other values fixed. But notice, if windspeed, radiation, pressure, and any feature in 
    the personal profiles change this lines and still this function could be used 
    in such cases ONLY if those variables are constant.
    
    
    Parameters
    ----------
    Survibability: boolean
        Matrix with the boolean results from the survivability assesment
    humidity : float
        Array with the humidities evaluated in the model 
    type_humidity : str
        Name of the humidity metric used in the analysis
    temperature : float
        Array with the humidities evaluate in the model
    tw : Pd.DataFrame
        Array with

        
    Returns
    -------
    New_survivability : pd.DataFrame'''  
                 
    New_survivability = pd.DataFrame(index = range(len(humidity)), columns = [type_humidity,'Tair','Tw'])
    New_survivability[type_humidity] = humidity
    
    #______________ Load data to the table ______________________________
    #Survivability
    for i,h_value in enumerate(humidity):
        index_array_h = np.where(humidity==h_value)[0][0]  
        survivability_row = Survivability[index_array_h,:]
        index_change_true_false = np.argwhere(np.diff(survivability_row)).squeeze()
        #DBT when survivability change from True to False given a rh level
        dbt_limit = temperature[index_change_true_false]    
        #Tw when survivability change from True to False given a rh level
        tw_limit = tw.iloc[index_array_h,index_change_true_false.tolist()]
        
        #Fill dataframe
        New_survivability['Tair'].loc[i] = dbt_limit
        New_survivability['Tw'].loc[i] = tw_limit
    
    return New_survivability


# =============================================================================
# Colormap 
# =============================================================================

import matplotlib.colors as clr

def Survivability_cmap():

    colores_cat =  [[252/255,252/255,205/255],
                    [237/255,232/255,131/255],
                    [254/255,192/255,169/255],
                    [205/255,206/255,254/255],
                    [191/255,200/255,209/255]] #Yellow, green, red, ,blue, grey
    
    cmap_survivability = clr.LinearSegmentedColormap.from_list('Surv', colores_cat, N=256)
    bins_survivability  = [0.5,1.5,2.5,3.5,4.5,5.5]
    norm_survivability = clr.BoundaryNorm(boundaries=bins_survivability, ncolors=256)
    ticks_survivability = [1,2,3,4,5]
    
    return cmap_survivability, bins_survivability, norm_survivability, ticks_survivability


    

