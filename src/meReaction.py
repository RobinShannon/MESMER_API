import numpy as np
from scipy.interpolate import BSpline
from scipy.interpolate import splrep
import matplotlib.pyplot as plt
from xml.dom import minidom
import math

class meReaction():

    def __init__(self, reacs, prods, ts, xml_root):
        self.reacs = reacs
        self.prods = prods
        self.ts = None
        self.ts = ts
        self.name = xml_root.attrib["id"]
        # TODO read other info from xml_root i.e tunneling and ILT details