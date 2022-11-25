import numpy as np
from scipy.interpolate import BSpline
from scipy.interpolate import splrep
import matplotlib.pyplot as plt
from lxml import etree as ET
import math

class meReaction():

    def __init__(self, reacs, prods, name, method='ILT', root=None, ts=None, ea = 0):
        self.reacs = reacs
        self.prods = prods
        self.ts = ts
        self.name = name
        self.tunneling = 'Eckart'
        self.method = method
        self.ilt_a = 1.0e-10
        self.ilt_n = 0
        self.ilt_ea = ea
        if root == None:
            self.cml = self.get_cml_snippet()
        else:
            self.cml = root

    def get_cml_snippet(self):
        ET.register_namespace('me', 'http://www.chem.leeds.ac.uk/mesmer')
        ET.register_namespace('xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        if len(self.prods) == 2 and len(self.reacs) == 1:
            temp = self.reacs
            self.reacs = self.prods
            self.prods = temp
        # Set up main me:mesmer xml tag and link to all the schemas
        reaction = ET.Element('reaction')
        reaction.set('id', str(self.name))
        reac = ET.SubElement(reaction, 'reactant')
        mol = ET.SubElement(reac, 'molecule')
        mol.set('role', 'modelled')
        mol.set('ref', str(self.reacs[0]))
        if len(self.reacs) == 2:
            reac2 = ET.SubElement(reaction, 'reactant')
            mol2 = ET.SubElement(reac2, 'molecule')
            mol2.set('role', 'excessReactant')
            mol2.set('ref', str(self.reacs[1]))
        prod = ET.SubElement(reaction, 'product')
        mol = ET.SubElement(prod, 'molecule')
        mol.set('role', 'modelled')
        mol.set('ref', str(self.prods[0]))

        if not self.ts == None:
            ts = ET.SubElement(reaction, 'transitionState')
            mol = ET.SubElement(ts, 'molecule')
            mol.set('role', 'transitionState')
            mol.set('ref', str(self.ts))
            method = ET.SubElement(reaction, '{http://www.chem.leeds.ac.uk/mesmer}MCRCMethod')
            method.set('{http://www.w3.org/2001/XMLSchema-instance}type','SimpleRRKM')
            tunneling =  ET.SubElement(reaction, '{http://www.chem.leeds.ac.uk/mesmer}tunneling')
            tunneling.set('{http://www.w3.org/2001/XMLSchema-instance}type','Eckart')
        else:
            method = ET.SubElement(reaction, '{http://www.chem.leeds.ac.uk/mesmer}MCRCMethod')
            method.set('{http://www.w3.org/2001/XMLSchema-instance}type','MesmerILT')
            A = ET.SubElement(method, '{http://www.chem.leeds.ac.uk/mesmer}preExponential')
            A.set('units',"cm3molecule-1s-1")
            A.text = '3e-10'
            n = ET.SubElement(method, '{http://www.chem.leeds.ac.uk/mesmer}nInfinity')
            n.text = '0'
            EA = ET.SubElement(method, '{http://www.chem.leeds.ac.uk/mesmer}activationEnergy')
            EA.set('units', "kJ/mol")
            EA.text = str(self.ilt_ea)
            excessReactant =  ET.SubElement(reaction, '{http://www.chem.leeds.ac.uk/mesmer}excessReactantConc')
            excessReactant.text = '1E16'
        return reaction






