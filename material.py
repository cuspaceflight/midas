import numpy as np

from tools import *

class Material(object):
    """Records a material for use in the propellant tanks"""
    def __init__(self, density, yield_stress):
        self.density = density
        self.yield_stress = yield_stress
    
    def calculate_tank_mass(self, pressure, radius, length, safety_factor):
        
        required_thickness = pressure * radius / (self.yield_stress / safety_factor)

        cylinder_mass = length * required_thickness * 2 * np.pi  * radius * self.density
        endcap_mass = 3 * np.pi * radius * required_thickness * self.density
        # Arbitraily assume the endcap is 3x the thickness, to account for fixing hardware plus extra structure. A hemisphere is effectively double

        mass = cylinder_mass + 2 * endcap_mass

        return(mass)