from tools import *

class Pressurant(object):
    """Object for handling a pressurant, initialised with just the type"""

    def __init__(self, pressurant_type):
        self.pressurant_type = pressurant_type

        if self.pressurant_type == Pressurant_Name.HELIUM:
            self.gamma = 1.67
        
        elif self.pressurant_type == Pressurant_Name.NITROGEN:
            self.gamma = 1.40
        
        elif self.pressurant_type == Pressurant_Name.CO2:
            self.gamma = 1.31

