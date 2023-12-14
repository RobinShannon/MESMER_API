import numpy as np
from scipy.interpolate import BSpline
from scipy.interpolate import splrep
import matplotlib.pyplot as plt
from xml.dom import minidom
import math

class meSpeciesProfile():

    def __init__(self, name, profile, T, P):
        self.name = name
        self.profile = profile
        self.T = T
        self.P = P

    def get_species(self, id):
        for i , mol in enumerate(self.name):
            if id == mol:
                return(self.profile[i])
