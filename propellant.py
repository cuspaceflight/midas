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
                self.Mr = 44

            properties = self.thermophys(self.temperature)
            self.den_l = properties["rho_l"]
            self.den_v = properties["rho_v"]
            self.h_l = properties["h_l"]
            self.h_v = properties["h_v"]
            self.pressure = properties["vapour_pressure"]

            self.fluid = Fluid(self.propellant_lookup["octopus_name"], method='helmholz', T = self.temperature, P = self.pressure)

            self.dryness = ((self.volume / self.mass) - (1/self.den_l)) / ((1/self.den_v) - (1/self.den_l))

            self.mass_l = (1 - self.dryness) * self.mass
            self.mass_v = self.dryness * self.mass
            self.enthalpy = self.mass_l * self.h_l + self.mass_v * self.h_v

            if self.mass_l/self.den_l > self.volume:
                print("Warning: too much mass of "+self.propellant_lookup["propep_name"])

        elif self.phase == Propellant_Phase.SELF_PRESSURISING and self.pressurant_name != None:
            self.pressurant = Pressurant(self.pressurant_name)
            self.target_pressure = pressure

            if self.propellant == Propellant_Name.NITROUS:
                self.thermophys = nitrous_thermophys
                self.Mr = 44
            
            if self.pressurant.pressurant_type == Pressurant_Name.HELIUM:
                self.pressurant_thermophys = helium_thermophys

            properties = self.thermophys(self.temperature)
            self.den_l = properties["rho_l"]
            self.den_v = properties["rho_v"]
            self.h_l = properties["h_l"]
            self.h_v = properties["h_v"]
            self.pressure = properties["vapour_pressure"]

            self.fluid = Fluid(self.propellant_lookup["octopus_name"], method='helmholz', T = self.temperature, P = self.pressure)

            self.dryness = ((self.volume / self.mass) - (1/self.den_l)) / ((1/self.den_v) - (1/self.den_l))

            self.mass_l = (1 - self.dryness) * self.mass
            self.mass_v = self.dryness * self.mass
            self.enthalpy = self.mass_l * self.h_l + self.mass_v * self.h_v

            if self.mass_l/self.den_l > self.volume:
                print("Warning: too much mass of "+self.propellant_lookup["propep_name"])
            
            pressure_deficit = self.target_pressure - self.pressure

            ullage_volume = self.volume - (self.mass_l/self.den_l)

            properties_p = self.pressurant_thermophys(self.temperature, self.target_pressure)
            self.den_p = properties_p["rho_v"]
            self.h_p = properties_p["h_v"]

            self.mass_p = ullage_volume * self.den_p * (pressure_deficit / self.target_pressure)

            self.pressure = self.target_pressure

            self.enthalpy += self.mass_p * self.h_p


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
                self.Mr = 60.096


            properties = self.thermophys(self.temperature)
            self.den_l = properties["rho_l"]
            self.den_v = properties["rho_v"]
            self.h_l = properties["h_l"]
            self.h_v = properties["h_v"]
            self.pressure = properties["vapour_pressure"]
            self.cv_l = properties["cv_l"]

            self.mass_l = prop_mass
            self.mass_v = 0

            if self.mass_l/self.den_l > self.volume:
                print("Warning: too much mass of "+self.propellant_lookup["propep_name"])

            if self.pressurant.pressurant_type == Pressurant_Name.HELIUM:
                self.pressurant_thermophys = helium_thermophys

            ullage_volume = self.volume - (self.mass_l/self.den_l)

            properties_p = self.pressurant_thermophys(self.temperature, self.target_pressure)
            self.den_p = properties_p["rho_v"]
            self.h_p = properties_p["h_v"]

            self.mass_p = ullage_volume*self.den_p

            self.enthalpy = self.mass_l*self.h_l + self.mass_p*self.h_p


        
        elif self.phase == Propellant_Phase.NONE:
            self.mass_l = self.mass_v = self.mass_p = 0.0

    def self_pressurising_equations_of_state_for_solving(self,arguments):
        """
        A wrapper onto the equations of state from thermophys that allows solving by scipy.optimise.fsolve
        arguments [float, float]:   [dryness fraction, temperature]
        """
        if arguments[1]>=309.57 or arguments[1]<=183.15 or arguments[0]<0 or arguments[0]>1:
            #If either parameter is out of bounds, return a random very large value to encourage the solver to move away
            H,V = 1e99,1e99
            return([H,V])
        else:
            dryness, temp = arguments
            properties = self.thermophys(temp)
            v_l = properties["v_l"]
            v_v = properties["v_v"]
            h_l = properties["h_l"]
            h_v = properties["h_v"]

            #enthalpy_diff and volume_diff are the difference between predicted and true values, squared to produce a smooth odd function
            if self.pressurant != None:
                enthalpy_diff = ((self.mass * dryness * h_v ) + (self.mass * (1 - dryness) * h_l) + (self.mass_p * self.h_p) - self.enthalpy)**2
                volume_diff = ((dryness * self.mass * v_v) + ((1-dryness) * self.mass * v_l) + (self.mass_p / self.den_p) - self.volume)**2
            else:
                enthalpy_diff = ((self.mass * dryness * h_v ) + (self.mass * (1 - dryness) * h_l) - self.enthalpy)**2
                volume_diff = ((dryness * self.mass * v_v) + ((1-dryness) * self.mass * v_l) - self.volume)**2

            return([enthalpy_diff,volume_diff])
    
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

    def remove_liquid(self,mass_out, add_pressurant_if_included=True):
        if mass_out > self.mass_l:
            return None

        if self.phase == Propellant_Phase.LIQUID and self.pressurant != None:   # Liquid state propellant
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

            properties_p = self.pressurant_thermophys(self.temperature, self.target_pressure)
            self.den_p = properties_p["rho_v"]
            self.h_p = properties_p["h_v"]

            self.enthalpy = self.mass_l*self.h_l + self.mass_p*self.h_p

        elif self.phase == Propellant_Phase.SELF_PRESSURISING: # Self-pressurising propellant with or without a gas supercharge added
            
            self.mass -= mass_out
            self.enthalpy -= mass_out * self.h_l
            
            #Estimate the temperature and dryness assuming no boiloff but homogenous temperature
            # The solver seems heavily biased towards the starting values so this should improve the overall results

            properties = self.thermophys(self.temperature)

            guess_temperature = (-self.h_l * mass_out) / ((properties["cp_l"]*self.mass_l) + (properties["cp_g"]*self.mass_v))
            guess_temperature += self.temperature
            if guess_temperature < 183.16:
                guess_temperature = 183.16

            guess_dryness = self.mass_v / (self.mass - mass_out)
            if guess_dryness>1:
                guess_dryness=1

            guess = [guess_dryness, guess_temperature]

            bounds = [(0,1),(183.16, 309.56)]

            out = least_squares(self.self_pressurising_equations_of_state_for_solving, guess, "3-point", bounds=bounds,xtol=1e-12,gtol=1e-12)

            self.dryness, self.temperature = out.x

            properties = self.thermophys(self.temperature)
            self.den_l = properties["rho_l"]
            self.den_v = properties["rho_v"]
            self.h_l = properties["h_l"]
            self.h_v = properties["h_v"]
            

            self.mass_l = (1 - self.dryness) * self.mass
            self.mass_v = self.dryness * self.mass

            if self.pressurant != None and add_pressurant_if_included:
                # Assume the pressurant is in thermal equilibrium with the rest of the propellant in the tank
                self.pressure = properties["vapour_pressure"]

                pressure_deficit = self.target_pressure - self.pressure

                ullage_volume = self.volume - (self.mass_l/self.den_l)

                properties_p = self.pressurant_thermophys(self.temperature, self.target_pressure)
                self.den_p = properties_p["rho_v"]
                self.h_p = properties_p["h_v"]

                required_mass_p = ullage_volume * self.den_p * (pressure_deficit / self.target_pressure)

                self.enthalpy += self.h_p * (required_mass_p - self.mass_p)
                self.mass_p = required_mass_p

                self.pressure = self.target_pressure

            elif self.pressurant != None:
                # No pressurant is added, the existing pressurant just expands
                # Note: to simplify the assumption, we're not tracking the temperature of the expanding pressurant
                # It's entirely constrained by the propellant temperature, mass and partial volume
                # This is fine during the run (where m_propellant ~ 1000 m_pressurant)
                # Towards the end of the run where the tank is mostly gas, it may be slightly inaccurate
                # Fix later I guess
                
                ullage_volume = self.volume - (self.mass_l/self.den_l)

                n_p = self.mass_p / self.pressurant.Mr
                n_v = self.mass_v / self.Mr

                total_pressure = properties["vapour_pressure"] * (n_p + n_v) / n_v

                self.pressure = total_pressure

            else:
                self.pressure = properties["vapour_pressure"]

    def heat_flux_through_wall(self, ambient_temperature, radius, dt, include_convection=False):
        # Quick calculation to account for heat flux through the wall of the propellant tank
        # when include_convection is true, this will probably become quite a hefty function
        # https://www.seas.upenn.edu/~lior/documents/HeattransferfromacylinderinaxialIJHMT.pdf
        # ^ Implement this 
        # Just a quick addition for now, to try to match Pulsar test data

        q = 0.5 * 5.67e-8 * (self.temperature**4 - ambient_temperature**4)

        area = 2 * self.volume / radius

        flux = q * area * dt

        self.enthalpy -= flux

        properties = self.thermophys(self.temperature)

        guess_temperature = (flux) / ((properties["cp_l"]*self.mass_l) + (properties["cp_g"]*self.mass_v))
        guess_temperature += self.temperature
        guess_dryness = self.dryness

        guess = [guess_dryness, guess_temperature]

        bounds = [(0,1),(183.16, 309.57)]

        out = least_squares(self.self_pressurising_equations_of_state_for_solving, guess, "3-point", bounds=bounds,xtol=1e-12,gtol=1e-12)

        self.dryness, self.temperature = out.x

        properties = self.thermophys(self.temperature)
        self.den_l = properties["rho_l"]
        self.den_v = properties["rho_v"]
        self.h_l = properties["h_l"]
        self.h_v = properties["h_v"]
        

        self.mass_l = (1 - self.dryness) * self.mass
        self.mass_v = self.dryness * self.mass





pulsar = Propellant("nitrous-self-pressurised", 270, 0.053, 33)
for i in range(0,18):
    pulsar.remove_liquid(1.4)
    pulsar.heat_flux_through_wall(275,0.1,1)
    print(pulsar.pressure/1e5, pulsar.temperature, pulsar.mass_l)