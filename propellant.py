from tools import *
import numpy as np
from scipy.optimize import *
from octopus.main import Fluid
from pressurant import Pressurant

class Propellant(object):
    """Object for handling a propellant, initialised without any idea of the properties 
    Arguments:
        temperature
        volume
        total mass
    
    Attributes:
        pressure
        dryness fraction 
        liquid density
        liquid enthalpy
        vapour density
        vapour enthalpy

    Functions:
        find_instrinsic_properties
            changing one of the arguments post-hoc, calculate all the others and then all attributes

    """

    def __init__(self, propellant, temperature, volume, prop_mass,  pressure = None):

        self.temperature = temperature
        self.volume = volume
        self.mass = prop_mass

        self.propellant_lookup = propellants[propellant]
        self.propellant = self.propellant_lookup["prop_name"]
        self.phase = self.propellant_lookup["phase"]
        self.pressurant_name = self.propellant_lookup["pressurant"]

        if self.phase == Propellant_Phase.SELF_PRESSURISING and self.pressurant_name == None:
            
            self.pressurant = None
            self.target_pressure = None
            self.mass_p = 0

            if self.propellant == Propellant_Name.NITROUS:
                self.thermophys = nitrous_thermophys

                nitrous_properties = self.thermophys(self.temperature)
                self.den_l = nitrous_properties["rho_l"]
                self.den_v = nitrous_properties["rho_v"]
                self.h_l = nitrous_properties["h_l"]
                self.h_v = nitrous_properties["h_v"]
                self.pressure = nitrous_properties["vapour_pressure"]

                self.fluid = Fluid('nitrous oxide', method='helmholz', T = self.temperature, P = self.pressure)

            self.dryness = ((self.volume / self.mass) - (1/self.den_l)) / ((1/self.den_v) - (1/self.den_l))

            

            self.mass_l = (1 - self.dryness) * self.mass
            self.mass_v = self.dryness * self.mass
            self.enthalpy = self.mass_l * self.h_l + self.mass_v * self.h_v

            if self.mass_l/self.den_l > self.volume:
                print("Warning: too much mass of "+self.propellant_lookup["propep_name"])

        elif self.phase == Propellant_Phase.LIQUID:
            if self.pressurant_name == None:
                raise RuntimeError
            else:
                self.pressurant = Pressurant(self.pressurant_name)

            if pressure != None:
                self.pressure = pressure
                self.target_pressure = pressure
            else:
                self.pressure = 10 * 1e5 

            if self.propellant == Propellant_Name.IPA:
                self.thermophys = ipa_thermophys

                ipa_properties = self.thermophys(self.temperature)
                self.den_l = ipa_properties["rho_l"]
                self.den_v = ipa_properties["rho_v"]
                self.h_l = ipa_properties["h_l"]
                self.h_v = ipa_properties["h_v"]
                self.pressure = ipa_properties["vapour_pressure"]
                self.cv_l = 2680

            self.mass_l = prop_mass
            self.mass_v = 0

            if self.mass_l/self.den_l > self.volume:
                print("Warning: too much mass of "+self.propellant_lookup["propep_name"])

            if self.pressurant.pressurant_type == Pressurant_Name.HELIUM:
                self.pressurant_thermophys = helium_thermophys
                self.den_p, self.h_p = self.pressurant_thermophys(self.temperature, self.pressure)

            self.mass_p = (self.volume - (self.mass_l/self.den_l))*self.den_p

            self.enthalpy = self.mass_l*self.h_l + self.mass_p*self.h_p
        
        elif self.phase == Propellant_Phase.NONE:
            self.mass_l = self.mass_v = self.mass_p = 0.0

    def self_pressurising_equations_of_state_for_solving(self,arguments,):
        """
        A wrapper onto the equations of state from thermophys that allows solving by scipy.optimise.fsolve
        arguments [float, float]:   [dryness fraction, temperature]
        """
        if arguments[1]>309.57 or arguments[1]<183.15 or arguments[0]<0 or arguments[0]>1:
            #If either parameter is out of bounds, return a random very large value to encourage the solver to move away
            H,V = 1e50,1e50
            return(H*V)
        else:
            dryness, temp = arguments
            properties = self.thermophys(self.temperature)
            v_l = properties["v_l"]
            v_v = properties["v_v"]
            h_l = properties["h_l"]
            h_v = properties["h_v"]

            #enthalpy_diff and volume_diff are the difference between predicted and true values, squared to produce a smooth odd function
            enthalpy_diff = ((self.mass * dryness * h_v ) + (self.mass * (1 - dryness) * h_l) - self.enthalpy)**2
            volume_diff = ((dryness * self.mass * v_v) + ((1-dryness) * self.mass * v_l) - self.volume)**2

            return((enthalpy_diff*volume_diff)**2)
    
    def liquid_pressurant_equations_of_state_for_solving(self,old_mass_l):
        """Provides a "solver" for the liquid-pressurant system

        Assumes the volume of IPA remains constant (note: this is not entirely valid but is handled in the iterative part of remove_liquid)

        Gas behaves adiabatically

        """
        gamma = self.pressurant.gamma


        old_volume_p = self.volume - (old_mass_l/self.den_l)
        new_volume_p = self.volume - (self.mass_l/self.den_l)

        new_temp = self.temperature * (old_volume_p/new_volume_p)**(gamma-1)

        new_pressure = self.pressure * (old_volume_p/new_volume_p)**(gamma) 

        return(new_temp, new_pressure)

    def remove_liquid(self, mass_out):

        if self.phase == Propellant_Phase.SELF_PRESSURISING and self.pressurant == None and self.propellant == Propellant_Name.NITROUS:
            self.mass -= mass_out
            self.enthalpy -= mass_out * self.h_l

            guess = np.array([self.dryness, self.temperature])

            bounds = [[0,1],[183.15,309.57]]

            new_dryness, new_temp = dual_annealing(self.self_pressurising_equations_of_state_for_solving,bounds,x0=guess, accept=-10.0, maxiter=1000).x

            self.dryness = new_dryness
            self.temperature = new_temp

            properties = self.thermophys(self.temperature)
            self.den_l = properties["rho_l"]
            self.den_v = properties["rho_v"]
            self.h_l = properties["h_l"]
            self.h_v = properties["h_v"]
            self.pressure = properties["vapour_pressure"]

            self.mass_l = (1 - self.dryness) * self.mass
            self.mass_v = self.dryness * self.mass
            self.enthalpy = self.mass_l * self.h_l + self.mass_v * self.h_v

        elif self.phase == Propellant_Phase.SELF_PRESSURISING and self.pressurant == None and self.propellant == Propellant_Name.NITROUS and False:
            self.mass -= mass_out
            self.enthalpy -= mass_out * self.h_l

            u = ['p', 'chi']
            dryness_estimate = self.dryness#self.mass_v / self.mass
            y = [self.pressure, self.dryness]

            initial = least_squares(self.fluid.fun_ps, [800, 250], bounds=([0, 0], [np.inf, self.fluid.Tc]), args=[u, y])

            den, new_temp = initial.x

            properties = self.fluid.get_properties(den, new_temp)

            self.dryness = properties["chi"]
            self.temperature = new_temp

            self.den_l, self.den_v, self.h_l, self.h_v, self.latent_heat_vapourisation, self.pressure, ldynvis = self.thermophys(self.temperature)
            self.mass_l = (1 - self.dryness) * self.mass
            self.mass_v = self.dryness * self.mass
            self.enthalpy = self.mass_l * self.h_l + self.mass_v * self.h_v


        elif self.phase == Propellant_Phase.LIQUID:

            self.mass_l -= mass_out
            self.enthalpy -= mass_out * self.h_l

            guess = [float(self.temperature)]

            bounds = [[180.0, 320.0]]

            for i in range(0,10):
                new_temp, new_pressure = self.liquid_pressurant_equations_of_state_for_solving(self.mass_l + mass_out)
                self.den_l = self.thermophys(new_temp)["rho_l"]

            self.temp, self.pressure = new_temp, new_pressure

            properties = self.thermophys(self.temperature)
            self.den_l = properties["rho_l"]
            self.den_v = properties["rho_v"]
            self.h_l = properties["h_l"]
            self.h_v = properties["h_v"]
            self.pressure = properties["vapour_pressure"]

            self.den_p, self.h_p = self.pressurant_thermophys(self.temperature, self.pressure)
            self.enthalpy = self.mass_l*self.h_l + self.mass_p*self.h_p

    def pressurant_input_equations_of_state_for_solving(self, temperature_in, pressure_in):
        """ Returns the added mass and new temperature when regulated pressurant is added to the system
        Source: CUED 1A Thermofluids Example paper 3, Question 6
        Bet you never thought you'd see an example paper cited in a real application!

        Warning: does not handle a lower pressure in regulator than tank. Currently assumes a check valve, ie will prevent mass flow out

        Arguments:
            temperature_in  float   input temperature of pressurant
            pressure_in     float   input pressure and target pressure of pressurant

        Returns:
            mass_in         float   mass of pressurant added
            final_temp      float   temperature of pressurant at end
        """
        gamma = self.pressurant.gamma
        coeffs = [gamma - 1 + (temperature_in/self.temperature),
                gamma - (pressure_in/self.pressure) + (temperature_in/self.temperature) + ((self.mass_l * self.cv_l) / (self.mass_p * 8314)),
                1 - (pressure_in/self.pressure) + ((self.mass_l * self.cv_l) / (self.mass_p * 8314)) * (1 - (pressure_in/self.pressure))]
        
        roots = np.roots(coeffs)
        if roots[0].imag !=0:
            roots = [0,0]
        
        f = max(roots)


        if f != 0:
            mass_in = f * self.mass_p

            final_temp = self.temperature * (pressure_in/self.pressure) / (1 + f)
        
        else:
            mass_in = 0
            final_temp = self.temperature

        return(mass_in, final_temp)

    def top_up_pressurant(self, input_temperature, target_pressure):
        """Uses an iterative solver

        """
        
        for i in range(0,10):
            mass_in, new_temp = self.pressurant_input_equations_of_state_for_solving(input_temperature, target_pressure)

        #new_temp, mass_in = dual_annealing(self.pressurant_input_equations_of_state_for_solving, bounds=[[180,320],[0,2]],args=[target_presure ,input_temperature], accept=-7.5, maxiter=2000).x

        self.temperature = new_temp
        self.mass_p += mass_in
        self.pressure = target_pressure

        properties = self.thermophys(self.temperature)
        self.den_l = properties["rho_l"]
        self.den_v = properties["rho_v"]
        self.h_l = properties["h_l"]
        self.h_v = properties["h_v"]
        self.pressure = properties["vapour_pressure"]

        self.den_p, self.h_p = self.pressurant_thermophys(self.temperature, self.pressure)
        self.enthalpy = self.mass_l*self.h_l + self.mass_p*self.h_p


nitrous = Propellant("nitrous-self-pressurised", 298, 0.1, 70)
nitrous.remove_liquid(0.1)