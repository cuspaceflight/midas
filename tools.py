import numpy as np
from enum import Enum


class Propellant_Phase(Enum):
    """Records the phase of material in a tank"""
    NONE = -1
    SELF_PRESSURISING       = 0
    LIQUID                  = 1

class Pressurant_Name(Enum):
    """Records the pressurant uses"""
    HELIUM      = 0
    NITROGEN    = 1
    CO2         = 2

class Propellant_Name(Enum):
    """Records the propellant used"""
    NONE = -1
    NITROUS = 0
    IPA = 10

propellants = {
    "nitrous-self-pressurised" : {
        "phase": Propellant_Phase.SELF_PRESSURISING,
        "prop_name": Propellant_Name.NITROUS,
        "pressurant": None,
        "propep_name": 'NITROUS OXIDE',
        "octopus_name": ' nitrous oxide'
    },
    "ipa-helium-pressurised": {
        "phase": Propellant_Phase.LIQUID,
        "prop_name": Propellant_Name.IPA,
        "pressurant": Pressurant_Name.HELIUM,
        "propep_name": 'ISOPROPYL ALCOHOL',
        "octopus_name": "isopropanol"
    },
    "none": {
        "phase": Propellant_Phase.NONE,
        "prop_name": Propellant_Name.NONE,
        "pressurant": None,
        "propep_name": 'NONE'
    }
}



def nitrous_thermophys(temp):
    """Get N2O data at a given temperature.
    Uses polynomials from ESDU sheet 91022. All units SI.
    Returns:
        - N2O liquid density
        - vapour density
        - liquid enthalphy
        - vapour enthalpy
        - latent heat of vaporization
        - dynamic viscosity
        - vapour pressure for input temperature
    """

    if not 183.15 <= temp <= 309.57:
        raise ValueError('nitrous oxide temperature out of data range')

    # Some handy definitions
    # ...I don't know what to call these properly
    T0 = temp / 309.57
    T0_RECIP = 1 / T0
    T0_INV = 1 - T0

    lden = 452 * np.exp(+ 1.72328 * pow(T0_INV, 1/3)
                        - 0.83950 * pow(T0_INV, 2/3)
                        + 0.51060 * T0_INV
                        - 0.10412 * pow(T0_INV, 4/3))

    vden = 452 * np.exp(- 1.009000 * pow(T0_RECIP - 1, 1/3)
                        - 6.287920 * pow(T0_RECIP - 1, 2/3)
                        + 7.503320 * (T0_RECIP - 1)
                        - 7.904630 * pow(T0_RECIP - 1, 4/3)
                        + 0.629427 * pow(T0_RECIP - 1, 5/3))

    hl = ((-200 + 116.043 * (T0_INV ** (1 / 3))+ -917.225 * (T0_INV ** (2 / 3))
           +794.779 * T0_INV + -589.587 * (T0_INV ** (4 / 3))) * 1000)

    hg = ((-200 + 440.055 * (T0_INV ** (1 / 3)) + -459.701 * ((1- (temp / 309.57)) ** (2 / 3))
           +434.081 * T0_INV + -485.338 * (T0_INV ** (4 / 3))) * 1000)

    c = ((2.49973 * (1 + 0.023454 * (T0_INV ** (-1)) + -3.80136 * T0_INV
                   + 13.0945 * (T0_INV ** 2) + -14.5180 * (T0_INV ** 3))) * 1000)

    vpres =  7251000 * (np.e **(T0_RECIP * 
        ( -6.71893 * T0_INV) + 
        (1.35966 * pow(T0_INV, 1.5)) + 
        (-1.3779 * pow(T0_INV, 2.5)) + 
        (-4.051 * pow(T0_INV, 5))))

    ldynvis = (0.0293423*np.e**((1.609*(((309.57-5.24)/(temp-5.24)-1)**(1/3)))
               +(2.0439*(((309.57-5.24)/(temp-5.24)-1)**(4/3)))))


    properties = {
        "rho_l" : lden,
        "rho_v" : vden,
        "h_l" :   hl,
        "h_v":    hg,
        "specific_heat_capacity_isobaric":   c,
        "vapour_pressure": vpres,
        "dynamic_visc": ldynvis,
        "v_l":  (1/lden),
        "v_v":  (1/vden)
    }
    return(properties)

def nitrous_thermophys_dt(temp):
    if not 183.15 <= temp <= 309.57:
        raise ValueError('nitrous oxide temperature out of data range')

    # Some handy definitions
    # ...I don't know what to call these properly
    T0 = temp / 309.57
    T_crit = 309.57
    T0_RECIP = 1 / T0
    T0_INV = 1 - T0
    
    nitrous_properties = nitrous_thermophys(temp)

    vpres_dT = nitrous_properties["vapour_pressure"] * (((T_crit*( 
            (6.71893 / T_crit) +
            ((-2.03949 * pow(T0_INV, 0.5)) / T_crit) + 
            (((3.44475 * pow(T0_INV, 1.5)) / T_crit) + 
            (20.255 * pow(T0_INV, 4)) / T_crit))
        / temp))  -
        ((T_crit / pow(temp, 2)) * (
            (-6.71893 * T0_INV) +
            (1.35966 * pow(T0_INV, 1.5)) +
            (-1.3779 * pow(T0_INV, 2.5)) + 
            (-4.051 * pow(T0_INV, 5)) 
            )))
    
    print(vpres_dT)

def ipa_thermophys(temp):
    """Get IPA data at a given temperature
    Uses polynomials from a range of sources. All units SI
    It would be really nice to get a hold of ESDU 89017 datasheet, alas he is expensive

    Format is designed to match nitrous_thermophys, hence the range of empty /redundant valyes

    Returns
        - IPA liquid density
        - vapour density (set at 0)
        - liquid enthalphy 
        - vapour enthalpy (set at 0)
        - latent heat of vaporization (set at 0)
        - dynamic viscosity
        - vapour pressure for input temperature

    """
    lden = 925.87 + (-0.0359 * temp) + (-0.0015 * temp**2)
    # Fit for data obtained from http://www.ddbst.com/en/EED/PCP/DEN_C95.php

    vden = 1e-7

    hl = -247.1 + (-0.6133 * temp) + (0.0055 * temp**2)
    # Fit for data obtained from https://pubs.acs.org/doi/10.1021/ie50466a033

    hg = 0

    c = 0

    vpres = 1000 * (20693 + (-119.57 * temp) + (0.1724 * temp**2))
    # Fit for data obtained from http://www.ddbst.com/en/EED/PCP/VAP_C95.php

    ldynvis = 22.397 * np.e ** (-0.031 * temp)
    # Fit for data obtained from http://www.ddbst.com/en/EED/PCP/VIS_C95.php

    properties = {
        "rho_l" : lden,
        "rho_v" : vden,
        "h_l" :   hl,
        "h_v":    hg,
        "latent_heat_vap":   c,
        "vapour_pressure": vpres,
        "dynamic_visc": ldynvis,
        "v_l":  (1/lden),
        "v_v":  (1/vden)
    }
    return(properties)


def helium_thermophys(temp, pressure):
    """Gets helium data at a given temperature and pressure with ideal gas assumptions

    Returns:
        - gas density
        - gas enthalpy
    """

    vden = pressure / (2080 * temp)

    hv = 5.19 * temp

    return (vden, hv)

nitrous_thermophys_dt(290.5)
