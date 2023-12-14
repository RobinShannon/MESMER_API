from src.Main import MESMER_API
import numpy as np
from scipy.optimize import least_squares as lm
from subprocess import Popen, PIPE

def get_chi( values, inputs, params):
    mol_dict = dict(zip(params,values))
    chi =[]
    for inp in inputs:
        inp.modify_in_place(mol_dict)
        #run mesmer
        p = Popen(['/Users/chmrsh/Documents/mesmerStoch/bin/mesmer','temp.xml'], stdout=PIPE, stderr=PIPE )
        stdout, stderr = p.communicate()
        out = stderr.decode("utf-8")
        out = MESMER_API()
        out.parse_me_xml('Mesmer_out.xml')
        chi += np.sum(out.get_chi_sq())
    return chi

def sumarise( values, inputs, params):
    mol_dict = dict(zip(params,values))
    results =[]
    for i,inp in enumerate(inputs):
        inp.modify_in_place(mol_dict)
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
params = ['TS_COC=O_CO[C]=O', 'TS_COC=O_[CH2]OC=O']
values = [-2.04236176,  1.57146613]


sol=lm(get_chi, values, method='trf',jac='3-point',max_nfev=2, args=(inputs,params))
Hess = np.linalg.inv(np.dot(sol.jac.T, sol.jac))
ave_error = (sol.fun).sum()/len(sol.fun)
dFit =np.sqrt(np.diag(Hess*ave_error))
print('variables = ' + str(sol.x))
print('chi values = ' + str(sol.fun))


for i,inp in enumerate(inputs):
    res = sumarise(values,[inp],params)
    with open(str(i)+'out', 'w') as fp:
        fp.write('Temp\texp\terr\tcal\tAld\tMeth\n')
        for item in res:
            fp.write(str(item[0])+'\t'+str(item[1])+'\t'+ str(item[2])+'\t'+ str(item[3])+'\t'+ str(item[4])+'\t'+ str(item[5])+"\n")