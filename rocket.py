from ambiance import Atmosphere

from scipy.optimize import *
from scipy import interpolate
import csv
import os.path

from tools import *
from engine import *
from propellant import *


class Rocket(object):
    def __init__(self, fuel, oxidiser, tank_material, engine, excess_dry_mass, radius):

        self.tank_material = tank_material

        if type(engine) == LiquidEngine:
            self.fuel = fuel
            self.oxidiser = oxidiser
            self.fuel_pressurant = fuel.pressurant
            self.oxidiser_pressurant = oxidiser.pressurant
            

            fuel_tank_length = self.fuel.volume / (np.pi * radius **2)
            ox_tank_length = self.oxidiser.volume / (np.pi * radius **2)
            """
            Safety factors:
                Fuel tank = 1.5
                Ox tank = 1.75
            """
            fuel_tank_SF = 1.5
            ox_tank_SF = 1.75

            self.fuel_tank_mass = self.tank_material.calculate_tank_mass(self.fuel.pressure, radius, fuel_tank_length, fuel_tank_SF)
            self.ox_tank_mass = self.tank_material.calculate_tank_mass(self.oxidiser.pressure, radius, ox_tank_length, ox_tank_SF)

        elif type(engine) == SolidMotor:
            self.fuel = Propellant("none", 300, 0, 0)
            self.oxidiser = Propellant("none", 300, 0, 0)

            self.fuel_tank_mass = self.ox_tank_mass = 0.0


        self.engine = engine
        self.excess_dry_mass = excess_dry_mass

        self.radius = radius

        self.mass = 0
        self.mass = self.fuel.mass_l + self.fuel.mass_v + self.fuel.mass_p + self.fuel_tank_mass 
        self.mass += self.oxidiser.mass_l + self.oxidiser.mass_v + self.oxidiser.mass_p + self.ox_tank_mass 
        self.mass += self.engine.mass + self.excess_dry_mass

        self.load_drag_data(1)
        self.load_drag_data(2)
    
    def thrust_timestep(self, dt=1, time=0.1, atm_pressure=1):

        if type(self.engine) == LiquidEngine:
            fuel_dP = self.fuel.pressure - self.engine.chamber_pressure
            ox_dP = self.oxidiser.pressure - self.engine.chamber_pressure

            if fuel_dP < 0 or ox_dP < 0:
                print("====== Main Engine Misfire")
                return(0)
            
            fuel_mass_flow, ox_mass_flow = self.engine.mass_flow_rate(self.oxidiser.pressure, self.oxidiser.temperature, self.fuel.pressure, self.fuel.temperature)

            if fuel_mass_flow * dt > self.fuel.mass_l or ox_mass_flow * dt > self.oxidiser.mass_l:
                print("===== Main Engine Cut Off =====")
                return(0)

            thrust = self.engine.get_thrust(fuel_mass_flow, ox_mass_flow, atm_pressure)

            
            self.fuel.remove_liquid(fuel_mass_flow * dt)
            self.oxidiser.remove_liquid(ox_mass_flow * dt)

            if self.fuel.pressurant != None: 
                self.fuel.top_up_pressurant(self.fuel.temperature, self.fuel.target_pressure)

            if self.oxidiser.pressurant != None: 
                self.fuel.top_up_pressurant(self.oxidiser.temperature, self.oxidiser.target_pressure)
        
        elif type(self.engine) == SolidMotor:
            thrust = self.engine.get_thrust(time)
            self.engine.get_mass(time)

        self.mass = self.fuel.mass_l + self.fuel.mass_v + self.fuel.mass_p + self.fuel_tank_mass 
        self.mass += self.oxidiser.mass_l + self.oxidiser.mass_v + self.oxidiser.mass_p + self.ox_tank_mass 
        self.mass += self.engine.mass + self.excess_dry_mass
        

        return(thrust)
    
    def load_drag_data(self, stage):
        """
        Reads the csv files with Cd-Mach data and returns them in a usable format.
        stage [int] should be 1 or 2 (because indexing starts from 1, obviously)

        Output is a pair of functions Cd_off(mach) and Cd_on(mach)
        """
        mach_data = []
        cd_on_data = []
        cd_off_data = []

        if stage == 1:
            filename = "CD_1st_stage"
        elif stage == 2:
            filename = "CD_2nd_stage"
        
        script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
        rel_path = filename+".CSV"
        abs_file_path = os.path.join(script_dir, rel_path)
        with open(abs_file_path) as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            linecount = 0
            for row in reader:
                if linecount!=0:
                    if row[1] == "0":
                        mach_data.append(float(row[0]))
                        cd_off_data.append(float(row[3]))
                        cd_on_data.append(float(row[4])) 
                linecount+=1


        self.cd_off = interpolate.interp1d(mach_data, cd_off_data)
        self.cd_on = interpolate.interp1d(mach_data, cd_on_data)
    
    def get_drag(self, altitude, speed, engine_on):
        if altitude >= 81020:
            return (0)
        atmosphere = Atmosphere(altitude)
        mach = speed/atmosphere.speed_of_sound[0]
        density = atmosphere.density[0]

        if mach <= 0.01:
            mach = 0.01
        if engine_on:
            cd = self.cd_on(mach)
        else:
            cd = self.cd_off(mach)
        
        

        drag_force = 0.5 * cd * density * (np.pi * self.radius ** 2) * (speed ** 2)

        return(drag_force)
