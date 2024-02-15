import numpy as np
from openbabel import openbabel, pybel
from ase import Atoms
from ase.io import write
from lxml import etree as ET
from io import StringIO
import sys


class meMolecule():

    def __init__(self, ase_mol, role = 'modelled', coupled=False, **kwargs):
        ET.register_namespace('me', 'http://www.chem.leeds.ac.uk/mesmer')
        ET.register_namespace('xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        self.ase_mol = ase_mol
        self.cml = self.write_geometry()
        self.name = kwargs.get('name', 'temp_name')
        self.zpe = kwargs.get('zpe', 0)
        self.newBonds = kwargs.get('newBonds', None)
        self.vibFreqs = kwargs.get('vibFreqs',[])
        self.hinderedRotors = kwargs.get('hinderedRotors', [])
        self.hinderedBonds = kwargs.get('hinderedBonds', [])
        self.hinderedAngles = kwargs.get('hinderedAngles', [])
        if self.newBonds != None:
            self.add_bonds(self.newBonds)
        atom_num = ase_mol.get_atomic_numbers()
        s = sum(atom_num)
        self.number_of_atoms = len(atom_num)
        if s % 2 == 0:
            default_spin  = 1
        else:
            default_spin = 2
        self.spinMultiplicity = kwargs.get('spinMultiplicity', default_spin)
        self.imaginary_freq = kwargs.get('imaginary_frequency', 0)
        self.hessian = kwargs.get('hessian', [])
        self.role = role
        self.epsilon = kwargs.get('epsilon', 500)
        self.sigma = kwargs.get('sigma', 5)
        self.de_down = kwargs.get('de_down', 200)
        try:
            self.add_properties()
        except:
            pass
        if coupled:
            self.add_hindered_rotor_coupled(self.hinderedRotors, self.hinderedAngles,self.hinderedBonds)
        else:
            for rot, bond, angle in zip(self.hinderedRotors, self.hinderedAngles,self.hinderedBonds):
                self.add_hindered_rotor(rot,bond,angle)

    @classmethod
    def from_cml(cls,cml, name, role = 'modelled'):
        data = {}
        data['name'] = name
        # using the xml element mol, get the details
        # first try to read the geometry
        try:
            ase_mol = cls.read_geometry(cml)
        except:
            print("couldnt read geometry for species " + str(name))
        # then loop through the properties
        props = cml.findall("{http://www.xml-cml.org/schema}propertyList")[0].findall("{http://www.xml-cml.org/schema}property")
        for prop in props:
            pid = prop.attrib["dictRef"]
            if pid == "me:ZPE":
                zpe = float(prop.findall("{http://www.xml-cml.org/schema}scalar")[0].text)
                data['zpe'] = zpe
            if pid == "me:vibFreqs":
                string = prop.findall("{http://www.xml-cml.org/schema}array")[0].text
                vib_freqs = [float(i) for i in string.split()]
                data['vib_freqs'] = vib_freqs
            if pid == "me:spinMultiplicity":
                spin = float(prop.findall("{http://www.xml-cml.org/schema}scalar")[0].text)
                data['spin'] = spin
            if pid == "me:epsilon":
                epsilon = float(prop.findall("{http://www.xml-cml.org/schema}scalar")[0].text)
                data['epsilon'] = epsilon
            if pid == "me:sigma":
                sigma = float(prop.findall("{http://www.xml-cml.org/schema}scalar")[0].text)
                data['sigma'] = sigma
        return cls(ase_mol, role, **data)



    @staticmethod
    def read_geometry(mol):
        # Create OBabel object from cml
        xml = ET.tostring(mol, encoding='unicode', method='xml')
        xml_head = '<?xml version="1.0" ?>\n<me:mesmer xmlns="http://www.xml-cml.org/schema" xmlns:me="http://www.chem.leeds.ac.uk/mesmer" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n'
        xml_bottom = '</me:mesmer>'
        xml = xml_head + xml + xml_bottom
        mol = pybel.readstring('cml', xml)
        dim = len(mol.atoms)
        a = np.zeros(dim)
        b = np.zeros((dim, 3))
        i = 0
        for Atom in mol:
            a[i] = Atom.atomicnum
            b[i] = Atom.coords
            i += 1
        ase_mol = Atoms(symbols=a, positions=b)
        return ase_mol

    def write_full_cml(self):
        geom = self.write_geometry()

    def add_bonds(self, newBonds):
        bonds = self.cml.findall("{http://www.xml-cml.org/schema}bondArray")[0]
        print(str(bonds))
        for changed in newBonds:
            b = ET.SubElement(bonds,'{http://www.xml-cml.org/schema}bond')
            b.set('atomRefs2','a' + str(int(changed[0])) + ' ' + 'a' +str(int(changed[1])))
            b.set('order', '1')
            b.text = None
        ET.indent(bonds, space='\n    ')

    def add_hindered_rotor(self, rotor_array, angle_array, bond_idx):
        # Set up main me:mesmer xml tag and link to all the schemas
        hind = ET.SubElement(self.cml,'{http://www.chem.leeds.ac.uk/mesmer}ExtraDOSCMethod')
        hind.set('{http://www.w3.org/2001/XMLSchema-instance}type', "me:HinderedRotorQM1D")

        bond_id = ET.SubElement(hind, '{http://www.chem.leeds.ac.uk/mesmer}bondRef')
        bond_id.text = 'bond'+str(bond_idx[0])+str(bond_idx[1])

        potential = ET.SubElement(hind, '{http://www.chem.leeds.ac.uk/mesmer}HinderedRotorPotential')
        hind.set('format', "numerical")
        hind.set('units', "kj/mol")
        hind.set('expansionSize', "10")
        inc = 360 / (len(rotor_array) -1)
        for h,a in zip(rotor_array, angle_array):
            point = ET.SubElement(potential, '{http://www.chem.leeds.ac.uk/mesmer}PotentialPoint')
            point.set('angle', str(a))
            point.set('potential', str(h))
        inert = ET.SubElement(hind, '{http://www.chem.leeds.ac.uk/mesmer}CalculateInternalRotorInertia')
        inert.set('phaseDifference', '0')

        # now find the appropriate bond in the list to mark

        bonds = self.cml.findall("{http://www.xml-cml.org/schema}bondArray")[0].findall("{http://www.xml-cml.org/schema}bond")
        for bond in bonds:
            pid = bond.attrib["atomRefs2"]
            if pid == str('a'+str(bond_idx[0]) + ' a' + str(bond_idx[1])) or pid == str('a'+str(bond_idx[1]) + ' a' + str(bond_idx[0])):
                bond.set('id', 'bond'+str(bond_idx[0])+ ' ' + str(bond_idx[1]) )
        ET.indent(hind, space = '\n    ')

    def add_hindered_rotor_coupled(self, rotor_array, angle_array, bond_idx):
        # Set up main me:mesmer xml tag and link to all the schemas
        couple = ET.SubElement(self.cml,'{http://www.chem.leeds.ac.uk/mesmer}ExtraDOSCMethod')
        couple.set('{http://www.w3.org/2001/XMLSchema-instance}type', '{http://www.chem.leeds.ac.uk/mesmer}ClassicalCoupledRotors')
        MC = ET.SubElement(couple,'{http://www.xml-cml.org/schema}MCPoints')
        MC.text='1000'
        RotArr = ET.SubElement(couple, '{http://www.chem.leeds.ac.uk/mesmer}RotorArray')
        for rot, angle, bond in zip(rotor_array, angle_array, bond_idx):
            hind = ET.SubElement(RotArr, '{http://www.chem.leeds.ac.uk/mesmer}Rotor')
            bond_string = 'bond'+str(bond[0])+ '_' +str(bond[1])
            hind.set('bondRef', bond_string)

            hind_pot = ET.SubElement(hind, '{http://www.chem.leeds.ac.uk/mesmer}HinderedRotorPotential')
            hind_pot.set('format', "numerical")
            hind_pot.set('units', "kJ/mol")
            hind_pot.set('expansionSize', "10")
            for h,a in zip(rot, angle):
                point = ET.SubElement(hind_pot, '{http://www.chem.leeds.ac.uk/mesmer}PotentialPoint')
                point.set('angle', str(a))
                point.set('potential', str(h))
            inert = ET.SubElement(hind, '{http://www.chem.leeds.ac.uk/mesmer}CalculateInternalRotorInertia')
            inert.set('phaseDifference', '0')

        # now find the appropriate bond in the list to mark
        bonds = self.cml.findall("{http://www.xml-cml.org/schema}bondArray")[0].findall("{http://www.xml-cml.org/schema}bond")
        for bond in bonds:
            pid = bond.attrib["atomRefs2"]
            if pid == str('a'+str(bond_idx[0]) + ' a' + str(bond_idx[1])) or pid == str('a'+str(bond_idx[1]) + ' a' + str(bond_idx[0])):
                bond.set('id', 'bond'+str(bond_idx[0])+ '_' + str(bond_idx[1]) )
        ET.indent(couple, space = '\n    ')

    def add_properties(self):
        # Set up main me:mesmer xml tag and link to all the schemas
        propList = ET.SubElement(self.cml,'propertyList')
        prop1 = ET.SubElement(propList, 'property')
        prop1.set('dictRef','me:ZPE')
        val1 = ET.SubElement(prop1,'scalar')
        val1.set('units', 'kJ/mol')
        val1.text = str(self.zpe)
        prop2 = ET.SubElement(propList, 'property')
        prop2.set('dictRef','me:vibFreqs')
        val2 = ET.SubElement(prop2,'array')
        val2.set('units', 'cm-1')
        val2.text = ' '.join(str(item) for item in self.vibFreqs)
        prop3 = ET.SubElement(propList, 'property')
        prop3.set('dictRef','me:spinMultiplicity')
        val3 = ET.SubElement(prop3,'scalar')
        val3.set('units', 'cm-1')
        val3.text = str(self.spinMultiplicity)
        prop5 = ET.SubElement(propList, 'property')
        prop5.set('dictRef','me:Hessian')
        val5 = ET.SubElement(prop5,'matrix')
        val5.set('matrixType', 'squareSymmetricLT')
        val5.set('rows', str(self.number_of_atoms*3))
        val5.set('units', 'kcal/mol/Ang2')
        val5.text = ' '.join(str(item) for item in self.hessian)
        if self.role == 'ts':
            prop4 = ET.SubElement(propList, 'property')
            prop4.set('dictRef', 'me:imFreqsScale')
            val4 = ET.SubElement(prop4, 'scalar')
            val4.set('units', 'cm-1')
            val4.text = str(self.imaginary_freq)
        ET.indent(propList, space = '\n    ')

    def write_geometry(self):
        # Create OBabel object from cml
        old_stdout = sys.stdout
        sys.stdout = xyz = StringIO()
        write('-', self.ase_mol, format='xyz')
        sys.stdout = old_stdout
        mol = pybel.readstring('xyz', xyz.getvalue())
        xyz.close()
        mol.addh()
        mol.make3D()
        str = mol.write(format = 'cml')
        lines = str.split('\n', 1)[-1]
        print(lines)
        cml = ET.fromstring(lines)
        return cml

    def write_cml(self, file = '-'):
        self.cml.set('id', str(self.name))
        ET.indent(self.cml, space='\n    ')
        xml = ET.tostring(self.cml, encoding='unicode', pretty_print=True)
        xml = xml.splitlines()
        strip_lines = [line for line in xml if line.strip() != ""]
        mod_xml= ""
        for line in strip_lines:
            mod_xml += line + "\n"
        f = open(file,"w")
        f.write(mod_xml)
        f.close()