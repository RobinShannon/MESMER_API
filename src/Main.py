import copy

from lxml import etree as ET
from src.meReaction import meReaction
from src.meMolecule import meMolecule
from src.mePT import mePT
from src.meSpeciesProfile import meSpeciesProfile


class MESMER_API():

    def __init__(self):
        # Initialise empty dicts for holding different sections of the xml
        self.reactions_dict = {}
        self.mols_dict = {}
        self.pt_list = []
        self.species_profile_list = []
        self.bartis_widom_dict = {}
        self.tree = ''

    def create_xml_shell(temperature = 298, pressure = 1):
        return

    def modify_in_place(self, energies_dict, imFreqs_dict=None):
        doc = self.tree.getroot()
        doc2 = copy.deepcopy(doc)
        mols = doc2.findall("{http://www.xml-cml.org/schema}moleculeList")[0].findall(
            "{http://www.xml-cml.org/schema}molecule")
        for mol in mols:
            nid = mol.attrib["id"]
            if nid in energies_dict.keys():
                zpe = energies_dict[nid]
                props = mol.findall("{http://www.xml-cml.org/schema}propertyList")[0].findall(
                    "{http://www.xml-cml.org/schema}property")
                for prop in props:
                    pid = prop.attrib["dictRef"]
                    if pid == "me:ZPE":
                        zpe_orriginal = float(prop.findall("{http://www.xml-cml.org/schema}scalar")[0].text)
                        zpe_modified = zpe_orriginal + zpe
                        prop.findall("{http://www.xml-cml.org/schema}scalar")[0].text = str(zpe_modified)
            if imFreqs_dict is not None and nid in imFreqs_dict.keys():
                imFreq = imFreqs_dict[nid]
                props = mol.findall("{http://www.xml-cml.org/schema}propertyList")[0].findall(
                    "{http://www.xml-cml.org/schema}property")
                for prop in props:
                    pid = prop.attrib["dictRef"]
                    if pid == "me:imaginaryFrequencyScaleFactor":
                        imFreq_orriginal = float(prop.findall("{http://www.xml-cml.org/schema}scalar")[0].text)
                        imFreq_modified = imFreq_orriginal + imFreq
                        prop.findall("{http://www.xml-cml.org/schema}scalar")[0].text = str(imFreq_modified)
        tree = ET.ElementTree(doc2)
        tree.write("temp.xml")

    def get_chi_sq(self):
        doc = self.tree.getroot()
        cond = doc.findall("{http://www.chem.leeds.ac.uk/mesmer}conditions")[0]
        PTs = cond.findall("{http://www.chem.leeds.ac.uk/mesmer}PTs")[0].findall("{http://www.chem.leeds.ac.uk/mesmer}PTpair")
        chi_total = []
        for pt in PTs:
            eig = pt.findall("{http://www.chem.leeds.ac.uk/mesmer}experimentalEigenvalue")
            err = eig[0].attrib["error"]
            calc = eig[0].attrib["calcVal"]
            exp = eig[0].text
            chi_total.append(abs(float(calc) - float(exp))/ float(err))
        return chi_total

    def get_eigen_vs_expt(self):
        doc = self.tree.getroot()
        cond = doc.findall("{http://www.chem.leeds.ac.uk/mesmer}conditions")[0]
        PTs = cond.findall("{http://www.chem.leeds.ac.uk/mesmer}PTs")[0].findall("{http://www.chem.leeds.ac.uk/mesmer}PTpair")
        eigs_array = cond.findall("{http://www.chem.leeds.ac.uk/mesmer}PTs")[0].findall("{http://www.chem.leeds.ac.uk/mesmer}experimentalEigenvalue")
        Ts = []
        exps = []
        errs = []
        calcs = []
        for pt in PTs:
            Ts.append(pt.attrib["T"])
            eigs= pt.findall("{http://www.chem.leeds.ac.uk/mesmer}experimentalEigenvalue")[0]
            errs.append(eigs.attrib["error"])
            calcs.append(eigs.attrib["calcVal"])
            exps.append(eigs.text)
        return [Ts,exps,errs,calcs]

    def parse_me_xml(self, iPath):
        # Parse mesmer xml to with minidom
        self.tree = ET.parse(iPath)
        doc = self.tree.getroot()
        try:
            self.read_species_profiles(doc)
        except:
            pass
        # Read the reactions first. Add reactions to to reactions_dict and then add the species from each reaction to the
        # mols dict
        self.read_reactions(doc, self.reactions_dict, self.mols_dict)

        # Then read the experimental conditions. Add the bath gas to the mols_dict.
        #self.read_pts(doc)
        # TODO Read the control section of the input
        # TODO Add any species not encountered in the reactions or as a bath gas to the mols_dict
        # If the xml is an output from a previous calculation then we can also read in results in the form of species
        # profiles or Bartis Widom rates
        #read_species_profiles(doc, species_profile_dict)
        #read_bartis_widom(doc, bartis_widom_dict)


    def read_reactions(self, doc, reactions_dict, mols_dict):
        # Navigate through the xml tree to the reactionList and get a list of reaction elements
        reactions = doc.findall("{http://www.xml-cml.org/schema}reactionList")[0].findall("{http://www.xml-cml.org/schema}reaction")

        # For each reaction get the reactants products and TS (if present)
        for reaction in reactions:
            # since there may be two reactants or products initialize empty lists to hold them
            r_names = []
            p_names = []
            reacs = reaction.findall("{http://www.xml-cml.org/schema}reactant")
            for r in reacs:
                # Get the name a role of each reactant, find the details of that reactant in the xml and add it to the
                # mols_dict
                name = (r.findall("{http://www.xml-cml.org/schema}molecule"))[0].attrib["ref"]
                role = (r.findall("{http://www.xml-cml.org/schema}molecule"))[0].attrib["role"]
                self.add_molecule_to_dict(doc, mols_dict, name, role)
                r_names.append(self.mols_dict[name])
            prods = reaction.findall("{http://www.xml-cml.org/schema}product")
            for p in prods:
                # Do the same for the products
                name = (p.findall("{http://www.xml-cml.org/schema}molecule"))[0].attrib["ref"]
                role = (p.findall("{http://www.xml-cml.org/schema}molecule"))[0].attrib["role"]
                self.add_molecule_to_dict(doc, mols_dict, name, role)
                p_names.append(self.mols_dict[name])
            t = reaction.findall("{http://www.chem.leeds.ac.uk/mesmer}transitionState")
            # And add the TS to the mols_dict if the reaction has one
            if len(t) > 0:
                name = (t[0].findall("{http://www.xml-cml.org/schema}molecule"))[0].attrib["ref"]
                self.add_molecule_to_dict(doc, mols_dict, name, 'ts')
                ts_name = self.mols_dict[name]
            else:
                ts_name = None
            # Instantiate a new reaction object and include the xml data to populate the remaining properties
            new_reac = meReaction(r_names, p_names, ts_name, reaction)
            # Add reaction object to the reactions_dict
            reactions_dict[new_reac.name] = new_reac


    def add_molecule_to_dict(self, doc, mols_dict, name, role):
        # Check if a species with the same name already exists in the mols_dict
        if name in mols_dict:
            return
        # If not go through the moleculeList checking for name
        mols = doc.findall("{http://www.xml-cml.org/schema}moleculeList")[0].findall("{http://www.xml-cml.org/schema}molecule")
        for mol in mols:
            nid = mol.attrib["id"]
            if nid == name:
                # Once the right element is found instantiate a new molecule object
                new_mol = meMolecule.from_cml(mol,name,role)
                # and add it to the dictionary
                mols_dict[name] = new_mol

    def read_pts(self, doc):
        pressure_temps = doc.findall("{http://www.chem.leeds.ac.uk/mesmer}conditions")[0].findall("{http://www.chem.leeds.ac.uk/mesmer}PTs")
        pts = pressure_temps[0].findall("{http://www.chem.leeds.ac.uk/mesmer}PTpair")
        for pt in pts:
            T = pt.attrib["{http://www.chem.leeds.ac.uk/mesmer}T"]
            P = pt.attrib["{http://www.chem.leeds.ac.uk/mesmer}T"]
            try:
                precision = pt.attrib["{http://www.chem.leeds.ac.uk/mesmer}T"]
            except:
                precision = 'd'
            newPT = mePT(T,P,precision)
            self.pt_list.append(newPT)

    def read_species_profiles(self, doc):
        species_profile = doc.findall("{http://www.chem.leeds.ac.uk/mesmer}analysis")[0].findall("{http://www.chem.leeds.ac.uk/mesmer}populationList")
        for sp in species_profile:
            T = sp.attrib["T"]
            P = sp.attrib["conc"]
            times = sp.findall("{http://www.chem.leeds.ac.uk/mesmer}population")
            species = times[0].findall("{http://www.chem.leeds.ac.uk/mesmer}pop")
            names = []
            profile = []
            for s in species:
                n = s.attrib["ref"]
                names.append(n)
                profile.append([])
            for t in times:
                species2 = t.findall("{http://www.chem.leeds.ac.uk/mesmer}pop")
                for i,s2 in enumerate(species2):
                    profile[i].append(float(s2.text))
            species_profile = meSpeciesProfile(names,profile,T,P)
            self.species_profile_list.append(species_profile)

    def read_BW_rates(self,doc):
        PTs = doc.findall("{http://www.chem.leeds.ac.uk/mesmer}PTs")[0].findall("{http://www.chem.leeds.ac.uk/mesmer}PTpair")
        for pt in PTs:
            eig = pt.findall("{http://www.chem.leeds.ac.uk/mesmer}experimentalEigenvalue")
            err = eig[0].attrib["error"]
            calc = eig[0].attrib["calculated"]
            exp = eig[0].text


