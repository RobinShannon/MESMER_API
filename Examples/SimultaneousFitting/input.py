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
        p = Popen(['/Users/chmrsh/Documents/MacMESMER/mesmer','temp.xml'], stdout=PIPE, stderr=PIPE )
        stdout, stderr = p.communicate()
        out = stderr.decode("utf-8")
        out = MESMER_API()
        out.parse_me_xml('Mesmer_out.xml')
        chi += out.get_chi_sq()
    return chi

me = MESMER_API()
me.parse_me_xml('NonDeut.xml')
me2 = MESMER_API()
me2.parse_me_xml('1D.xml')
out = MESMER_API()
inputs = [me,me2]
params = ['TS_COC=O_CO[C]=O', 'TS_COC=O_[CH2]OC=O']
values = [0,0]

sol=lm(get_chi, values, method='lm',diff_step=0.1, args=(inputs,params))
print('all done')