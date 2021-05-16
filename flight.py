import csv
import os.path
from scipy import interpolate
import ambiance
import matplotlib.pyplot as plt


from engine import *
from propellant import Propellant
from pressurant import Pressurant
from tools import *
from material import Material
from rocket import Rocket

def make_panthera(temperature, of_ratio , ox_tank_volume = None, ox_mass = 110):
    # A test with 30s burn time on White Giant
    
    aluminium = Material(2700, 270e6)

    nitrous_density = nitrous_thermophys(temperature)["rho_l"]
    ipa_density = ipa_thermophys(temperature)["rho_l"]

    # Densities multiplied by 0.9 to allow for ullage space
    if ox_tank_volume == None and (type(ox_mass) == int or type(ox_mass) == float):
        # We are doing a mass-driven rocket. Mass of propellant is fixed, size is not
        ox_tank_volume = ox_mass / (nitrous_density * 0.9)

        fuel_mass = ox_mass * (1 / (1 + of_ratio))
        fuel_tank_volume = fuel_mass / (ipa_density * 0.9)
    
    elif ox_mass == None and (type(ox_tank_volume) == int or type(ox_tank_volume) == float):
        # We are doing a volume-driven rocket. Volume of tanks (and thus) height of tanks) is fixed, mass is not
        fuel_mass = ox_mass * (1 / (1 + of_ratio))

        ox_mass = ox_tank_volume * (nitrous_density * 0.9)
        fuel_tank_volume = fuel_mass / (ipa_density * 0.9)

    
    nitrous = Propellant("nitrous-self-pressurised", temperature, ox_tank_volume, ox_mass)
    ipa = Propellant("ipa-helium-pressurised", temperature, fuel_tank_volume, fuel_mass, nitrous.pressure)

    # Recommended to have a 20-25% pressure drop across the injector, rounded up to 30% to account for additional plumbing pressure losses
    # This is equal to a loss of 25% of upstream pressure
    # Replace this with actual design data at a later date

    chamber_pressure = nitrous.pressure * 0.5

    target_ox_flow = 3.25
    target_fuel_flow = target_ox_flow * (1 / (1 + of_ratio))

    white_giant = LiquidEngine("ipa-helium-pressurised", "nitrous-self-pressurised", chamber_pressure, target_ox_flow+target_fuel_flow, of_ratio, 4, 0.15)
    #white_giant.construct_injector(nitrous.pressure, nitrous.temperature, ipa.pressure, ipa.temperature)


    panthera = Rocket(ipa, nitrous, aluminium, white_giant, 10, 0.15)
    return(panthera)


def make_two_stage_rocket(temperature, ox_mass, ox_volume = None):
    panthera = make_panthera(temperature, 3.5, ox_tank_volume=None, ox_mass=ox_mass)
    panthera.load_drag_data(1)
    #panthera.engine.construct_injector(panthera.oxidiser.pressure, panthera.oxidiser.temperature, panthera.fuel.pressure, panthera.fuel.temperature)

    pro98 = SolidMotor("Pro98", "white-dwarf")
    pro98.read_eng_file()
    aluminium = Material(2700, 270e6)

    condor = Rocket(None, None, aluminium, pro98, 15, 0.1)
    condor.load_drag_data(2)

    rocket = {
        "booster": panthera,
        "booster-sep" : False,
        "delay": 1.0,
        "sustainer": condor
    }

    return(rocket)

def integrator(x, v, a, dt):
    """
    Semi-implicit Euler method

    Leave me alone, it works. If you want RK, do it yourself
    """
    new_v = v + (a * dt)
    new_x = x + (new_v * dt)

    return(new_x, new_v)

def main_flight_loop(rocket, dt = 0.1):
    time = 0

    times = []

    position_1st_stage = []
    velocity_1st_stage = []
    acceleration_1st_stage = []

    position_2nd_stage = []
    velocity_2nd_stage = []
    acceleration_2nd_stage = []

    x_1st = 0
    v_1st = 0

    x_2nd = x_1st
    v_2nd = v_1st

    acc_1st = acc_2nd = 0

    stage1_burnout = 9999
    stage2_ignition = 9999
    stage2_burnout = 9999

    for i in range(int(200/dt)):

        times.append(time)

        position_1st_stage.append(x_1st)
        velocity_1st_stage.append(v_1st)
        acceleration_1st_stage.append(acc_1st)

        position_2nd_stage.append(x_2nd)
        velocity_2nd_stage.append(v_2nd)
        acceleration_2nd_stage.append(acc_2nd)
        
        if x_1st <= 81020: 
            atm_1st = ambiance.Atmosphere(x_1st)
            pressure_1st = atm_1st.pressure[0] / 1e5
        else: pressure_1st = 0

        if x_2nd <= 81020:
            atm_2nd = ambiance.Atmosphere(x_2nd)
            pressure_2nd = atm_2nd.pressure[0] / 1e5
        else: pressure_2nd = 0

        if not rocket["booster-sep"]:
            # If the booster is still attached
            thrust = rocket["booster"].thrust_timestep(dt, time, pressure_1st)

            drag = rocket["booster"].get_drag(x_1st, abs(v_1st), thrust!=0) * np.sign(v_1st)

            mass = rocket["booster"].mass + rocket["sustainer"].mass

            acceleration = (thrust - drag - mass * 9.81) / mass
            acc_1st = acceleration
            acc_2nd = acceleration

            #print(time, thrust, mass, acceleration)

            x_1st, v_1st = integrator(x_1st, v_1st, acceleration, dt)

            x_2nd = x_1st
            v_2nd = v_1st
            
            if thrust == 0 and stage1_burnout == 9999:
                stage1_burnout = time
        
        elif rocket["booster-sep"]:
            # If the sustainer has fired

            # Sustainer loop
            if x_2nd > 0:
                thrust = rocket["sustainer"].thrust_timestep(dt, time - stage2_ignition, pressure_2nd)

                drag = rocket["sustainer"].get_drag(x_2nd, abs(v_2nd), thrust!=0) * np.sign(v_2nd)

                mass = rocket["sustainer"].mass

                acceleration = (thrust - drag - mass * 9.81) / mass

                acc_2nd = acceleration

                x_2nd, v_2nd = integrator(x_2nd, v_2nd, acceleration, dt)
            else:
                print("a")
                break

            # Spent booster loop
            if x_1st > 0:
                drag = rocket["booster"].get_drag(x_1st, abs(v_1st), False) * np.sign(v_1st)

                mass = rocket["booster"].mass

                acceleration = ( - drag - mass * 9.81) / mass

                acc_1st = acceleration

                x_1st, v_1st = integrator(x_1st, v_1st, acceleration, dt)

        
        if stage1_burnout + rocket["delay"] <= time and not rocket["booster-sep"]:
            # After booster stage burnout occurs, separate the booster and fire the sustainer stage
            rocket["booster-sep"] = True
            stage2_ignition = time
        
        time += dt

        #print(time, x_1st, x_2nd)


    return(times, position_1st_stage, velocity_1st_stage, position_2nd_stage, velocity_2nd_stage, acceleration_1st_stage, acceleration_2nd_stage)

def plot_stages(times, position_1st_stage, velocity_1st_stage, position_2nd_stage, velocity_2nd_stage, acceleration_1st_stage, acceleration_2nd_stage):


    print("Maximum 2nd stage altitude: {:.2f}".format(max(position_2nd_stage)))
    
    fig, axs = plt.subplots(3, 1)

    axs[0].plot(times, position_1st_stage, times, position_2nd_stage)

    axs[1].plot(times, velocity_1st_stage, times, velocity_2nd_stage)

    axs[2].plot(times, acceleration_1st_stage, times, acceleration_2nd_stage)

    plt.show()


## This is currently configured for the timestep consistency experiment. It takes ages to run. Don't!

"""
t_steps = [1.0,0.75,0.5,0.25]
heights = [[],[],[],[]]

for i in range(len(t_steps)):
    t_step = t_steps[i]
    for j in range(0,10):
        griffin = make_two_stage_rocket(250,110)
        times, position_1st_stage, velocity_1st_stage, position_2nd_stage, velocity_2nd_stage, acceleration_1st_stage, acceleration_2nd_stage = main_flight_loop(griffin, t_step)
        max_alt = np.round(max(position_2nd_stage),1)
        heights[i].append(max_alt)
        print(t_step, max_alt)
        #plot_stages(times, position_1st_stage, velocity_1st_stage, position_2nd_stage, velocity_2nd_stage, acceleration_1st_stage, acceleration_2nd_stage)

print(heights)
"""
griffin = make_two_stage_rocket(260,110)

print(griffin["booster"].oxidiser.pressure/1e5)
print(griffin["booster"].mass)
print(griffin["booster"].oxidiser.mass)
print(griffin["booster"].fuel.mass)

#times, position_1st_stage, velocity_1st_stage, position_2nd_stage, velocity_2nd_stage, acceleration_1st_stage, acceleration_2nd_stage = main_flight_loop(griffin, 1)
        