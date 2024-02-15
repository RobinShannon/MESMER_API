from src.Main import MESMER_API
import numpy as np
from lmfit import Model, Parameters
from subprocess import Popen, PIPE

def get_chi( inputs, params):
    chi =[]
    for inp in inputs:
        inp.modify_in_place(params)
        #run mesmer
        p = Popen(['/Users/chmrsh/Documents/mesmerStoch/bin/mesmer','temp.xml'], stdout=PIPE, stderr=PIPE )
        stdout, stderr = p.communicate()
        out = stderr.decode("utf-8")
        out = MESMER_API()
        out.parse_me_xml('Mesmer_out.xml')
        chi.append(out.get_chi_sq())
    return chi


def sumarise( inputs, params):
    results =[]
    for i,inp in enumerate(inputs):
        inp.modify_in_place(params)
        #run mesmer
        p = Popen(['/Users/chmrsh/Documents/mesmerStoch/bin/mesmer','temp.xml'], stdout=PIPE, stderr=PIPE )
        stdout, stderr = p.communicate()
        out = stderr.decode("utf-8")
        out = MESMER_API()
        out.parse_me_xml('mesmer_out.xml')
        eigs = out.get_eigen_vs_expt()
        for i,sp in enumerate(out.species_profile_list):
            list1 = sp.get_species('CH_post_comp')
            list2 = sp.get_species('CH3_post_comp')
            branch1 = list1[-1]
            branch2 = list2[-1]
            results.append([eigs[0][i],eigs[1][i],eigs[2][i],eigs[3][i],branch1,branch2])
    return results

me = MESMER_API()
me.parse_me_xml('NonDeut.xml')
me2 = MESMER_API()
me2.parse_me_xml('1D.xml')
me3 = MESMER_API()
me3.parse_me_xml('3D.xml')
me4 = MESMER_API()
me4.parse_me_xml('4D.xml')
me5 = MESMER_API()
me5.parse_me_xml('OD.xml')
me6 = MESMER_API()
me6.parse_me_xml('1DOD.xml')

inputs = [me,me2]
values = [-2.04236176,  1.57146613]

# Create an lmfit Model
model = Model(get_chi)

# Set initial parameter values
params = Parameters()
params.add('TS_COC=O_CO[C]=O', value=-2.04236176)
params.add('TS_COC=O_[CH2]OC=O', value=1.57146613)

# Perform the Levenberg-Marquardt minimization
result = model.fit([0,0], params, x=inputs)

# Print the fit results
print(result.fit_report())

# Print uncertainties on the fitted values
uncertainties = np.sqrt(np.diag(result.covar))
for name, value, uncertainty in zip(result.var_names, result.best_values.values(), uncertainties):
    print(f'{name}: {value} ± {uncertainty}')