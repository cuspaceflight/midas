import pypropep as ppp 
from tools import *
from scipy.optimize import *
from scipy import interpolate
import csv
import os.path
from ambiance import Atmosphere

from octopus import Orifice, Manifold, Fluid, PropertySource

ppp.init()


class LiquidEngine(object):
    def __init__(self, fuel, ox, nominal_chamber_pressure, nominal_mass_flow, nominal_of, expansion_ratio, nozzle_radius):
        self.fuel_reference = propellants[fuel]
        self.ox_reference = propellants[ox]

        self.fuel = ppp.PROPELLANTS[self.fuel_reference["propep_name"]]
        self.ox = ppp.PROPELLANTS[self.ox_reference["propep_name"]]

        self.nominal_chamber_pressure = nominal_chamber_pressure
        self.chamber_pressure = nominal_chamber_pressure
        self.nominal_mass_flow = nominal_mass_flow
        self.nominal_of = nominal_of

        self.expansion_ratio = expansion_ratio
        self.nozzle_radius = nozzle_radius

        self.n_ox_orifices = 1
        self.n_fuel_orifices = 1

        self.performance = ppp.ShiftingPerformance()

        self.calculate_dry_mass()

    def thrust_nominal(self):
        if self.fuel_reference["propep_name"] == 'ISOPROPYL ALCOHOL' and self.ox_reference["propep_name"] == 'NITROUS OXIDE':
            OF = 3.5 

        return(self.get_thrust(self.nominal_mass_flow * (1 / (OF + 1)), self.nominal_mass_flow * (OF / (OF + 1)), 1))

    def get_thrust(self, fuel_mass_flow, ox_mass_flow, atm_pressure):
        """
        Assumes a straight linear scaling between chamber pressure and mass flow
        """
        mass_flow = fuel_mass_flow + ox_mass_flow
        rated_performance = mass_flow/self.nominal_mass_flow
        self.chamber_pressure = self.nominal_chamber_pressure * rated_performance

        self.performance = ppp.ShiftingPerformance()

        self.performance.add_propellants_by_mass([(self.fuel, fuel_mass_flow),(self.ox, ox_mass_flow)])

        self.performance.set_state(P = self.chamber_pressure, Ae_At=self.expansion_ratio)

        exhaust_velocity = self.performance.performance.Isp
        exhaust_pressure = self.performance.properties[2].P 
        c_star = self.performance.performance.cstar

        pressure_adjusted_isp = exhaust_velocity + (c_star * self.expansion_ratio) * ((exhaust_pressure / self.chamber_pressure) - (atm_pressure / self.chamber_pressure))

        thrust = pressure_adjusted_isp * mass_flow

        return(thrust)

    def calculate_dry_mass(self):
        # This is just a dark magic estimate. Update with some realistic numbers at some point
        self.mass = 100 * self.nozzle_radius**3 * (self.nominal_chamber_pressure/1e5) + 1
        
        return(self.mass)
    
    def construct_injector(self, oxidiser_pressure, oxidiser_temperature, fuel_pressure, fuel_temperature):
        target_fuel_flow = self.nominal_mass_flow * (1 / (1 + self.nominal_of))
        target_ox_flow = self.nominal_mass_flow * (self.nominal_of / (1 + self.nominal_of))

        fuel_flow, ox_flow = self.mass_flow_rate(oxidiser_pressure, oxidiser_temperature, fuel_pressure, fuel_temperature)

        self.n_fuel_orifices = int(target_fuel_flow / fuel_flow)
        self.n_ox_orifices = int(target_ox_flow / ox_flow)

    def mass_flow_rate(self, oxidiser_pressure, oxidiser_temperature, fuel_pressure, fuel_temperature):
        """
        Uses the beautiful Octopus to get a mass flow rate using actual injector equations.

        Dimensions of manifold etc are based on Octopus example case

        TODO: don't decide what equation sets are used for the ox/fuel. Actually check somewhere (idk where?)
        """

        oxidiser = Fluid(self.ox_reference["octopus_name"], P=oxidiser_pressure, method="helmholz")
        fuel = Fluid(self.fuel_reference["octopus_name"], P=fuel_pressure, T=fuel_temperature, method="thermo")

        ox_manifold = Manifold(oxidiser, PropertySource(p=oxidiser_pressure, T=oxidiser_temperature))
        fuel_manifold = Manifold(fuel, PropertySource(p=fuel_pressure, T=fuel_temperature))

        ox_orifice = Orifice(ox_manifold, 1e-2, 2e-3, orifice_type=0)
        fuel_orifice = Orifice(fuel_manifold, 1e-2, 1e-3, orifice_type=0)

        print(self.chamber_pressure)
        fuel_mass_flow = fuel_orifice.m_dot_SPI(self.chamber_pressure) * self.n_fuel_orifices
        dyer_flow_rate = ox_orifice.m_dot_dyer(self.chamber_pressure)
        if dyer_flow_rate!=None:
            ox_mass_flow = dyer_flow_rate * self.n_ox_orifices
        else:
            print(self.chamber_pressure, oxidiser_pressure)
            ox_mass_flow = self.nominal_mass_flow * (self.nominal_of)/(1 + self.nominal_of)

        return(fuel_mass_flow, ox_mass_flow)

class SolidMotor(object):
    def __init__(self, name, eng_file):
        self.name = name
        self.eng_file = eng_file
    
    def read_eng_file(self):
        """
        Takes a .eng motor performance file (the standard for engine data) and processes it
        Returns self.thrust_curve and self.mass_curve which are both functions (generated by scipy.interpolate)
        which take a time (since ignition) and return the thrust and mass of the motor at that time

        Data on .eng formats comes from https://www.thrustcurve.org/info/raspformat.html
        """

        script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
        rel_path = self.eng_file+".eng"
        abs_file_path = os.path.join(script_dir, rel_path)

        with open(abs_file_path) as eng_data:
            entries = eng_data.read().split(";\n")
        split_entries = []
        for i in entries:
            #Filter out blank entries, which usually occur when multiple returns are put together
            if len(i) > 5:
                split_entries.append(i.split("\n"))

        for i in split_entries:
            if self.name in i[0]:
                raw_eng_data = []
                for j in range(len(i)):

                    #This removes any empty entries from the string - for example if there's a string of spaces
                    split = [x for x in i[j].split(" ") if (x and x!=";")]
                    if len(split) != 0:
                        raw_eng_data.append(split)
        
        # Identify the header line of the motor entry
        header_line = -1
        start_of_data = -1

        for i in range(len(raw_eng_data)):
            j = raw_eng_data[i]
            if (len(j) == 7):
                # Make absolutely sure we haven't just a comment with 7 words. If you have a comment with 7 words that fits
                # this pattern, god help you! 
                if ((not j[0].isnumeric()) and (j[1].isnumeric()) and (j[2].isnumeric()) or ((not j[3].isnumeric()) 
            or (j[3] == "0")) and (j[4].isnumeric()) and (j[5].isnumeric()) and (j[6].isnumeric())):
                    header_line = i
                    start_of_data = i + 1


        # Pull the relevant data out of the header
        self.nozzle_radius = float(raw_eng_data[header_line][1]) / 2
        self.length = float(raw_eng_data[header_line][2])

        self.wet_mass = float(raw_eng_data[header_line][5])
        self.dry_mass = self.wet_mass - float(raw_eng_data[header_line][4])

        self.mass = self.wet_mass

        prop_mass = self.wet_mass - self.dry_mass

        self.burn_time = float(raw_eng_data[-1][0])

        # Pull the time and thrust curve out of the rest of the data
        time_data = []
        thrust_data = []
        mass_data = []

        for i in raw_eng_data[start_of_data:]:
            time_data.append(float(i[0]))
            thrust_data.append(float(i[1]))
        
        self.impulse = np.trapz(thrust_data, time_data)

        # Assuming constant Isp, find the mass of the motor at each stage of the burn
        for i in range(len(raw_eng_data[start_of_data:])):
            impulse_so_far = np.trapz(thrust_data[0:i], time_data[0:i])
            live_prop_mass = (1 - (impulse_so_far/self.impulse)) * prop_mass
            mass_data.append(self.dry_mass + live_prop_mass)
        
        self.thrust_curve_interp = interpolate.interp1d(time_data, thrust_data)
        self.mass_curve_interp = interpolate.interp1d(time_data, mass_data)

    def get_thrust(self, time):
        if time <= 0 or time >= self.burn_time:
            return(0)
        else:
            return(self.thrust_curve_interp(time))
    
    def get_mass(self, time):
        if time <= 0:
            self.mass = self.wet_mass
        elif time >= self.burn_time:
            self.mass = self.dry_mass
        else:
            self.mass = self.mass_curve_interp(time)
        




#panthera = make_panthera(295, 3.5)
