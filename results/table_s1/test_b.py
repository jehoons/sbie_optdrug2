import json,itertools,pytest
from os.path import exists
from numpy.random import random
from boolean3_addon import attr_cy, to_logic
from pdb import set_trace
from boolean3_addon import to_logic 
from os.path import exists,dirname
import pandas as pd
from tqdm import tqdm 

from sbie_optdrug.result_season2 import tab_s1
from sbie_optdrug.result_season2.tab_s1 import fumia_phenotype
from fumia_phenotype import *

# @pytest.mark.skip(reason='')
# def test_1():
#     import engine 
#     results = []    
#     input_string = '00000'
#     for i,s in enumerate(['State_Mutagen', 'State_GFs','State_Nutrients',
#                             'State_TNFalpha','State_Hypoxia']):    
#         template[s] = True if input_string[i] == '1' else False

#     inputs = [ template ]

#     for sim_input in inputs:
#         on_states = [] 
#         off_states = [] 
#         for lbl in sim_input:
#             if sim_input[lbl] == True: 
#                 on_states.append(lbl)
#             else: 
#                 off_states.append(lbl)

#         res = engine.main(samples=1000, steps=60, debug=False, 
#                             on_states=on_states, off_states=off_states)

#         res = attach_phenotype({'input_condition':sim_input, 'simul_result': res})
#         results.append(res)

#     json.dump(results, open(file_test1, 'w'), indent=4)


# @pytest.mark.skip(reason='')
# def test_2():
#     import engine 

#     results = []    
#     input_string = '11100'
#     for i,s in enumerate(['State_Mutagen', 'State_GFs','State_Nutrients',
#                             'State_TNFalpha','State_Hypoxia']):    
#         template[s] = True if input_string[i] == '1' else False

#     inputs = [ template ]

#     for sim_input in inputs:
#         on_states = [] 
#         off_states = [] 
#         for lbl in sim_input:
#             if sim_input[lbl] == True: 
#                 on_states.append(lbl)
#             else: 
#                 off_states.append(lbl)
#         # res = attr_cy.run(samples=1000, steps=60, debug=False, 
#         #                     on_states=on_states, off_states=off_states)
#         res = engine.main(steps=60, samples=1000, debug=False, 
#                 on_states=on_states, off_states=off_states)

#         res = attach_phenotype({'input_condition':sim_input, 'simul_result': res})
#         results.append(res)

#     json.dump(results, open(file_test2, 'w'), indent=4)


# @pytest.mark.skip(reason='')
def test_all(force, progress):
    if exists(file_b1) and force == False:
        return;

    import engine 

    inputs = [ {
        'State_Mutagen' : State_Mutagen, 
        'State_GFs': State_GFs,
        'State_Nutrients': State_Nutrients,
        'State_TNFalpha': State_TNFalpha,
        'State_Hypoxia': State_Hypoxia,
        'State_Gli': False,
        } for State_Mutagen,State_GFs,State_Nutrients,State_TNFalpha, \
            State_Hypoxia in itertools.product([False,True], repeat=5)]

    results = []
    steps = 60
    samples = 1000
    repeats = 100
    for sim_input in tqdm(inputs, ascii=True):
        on_states = [] 
        off_states = [] 
        for lbl in sim_input:
            if sim_input[lbl] == True: 
                on_states.append(lbl)
            else: 
                off_states.append(lbl)

        res = attr_cy.parallel(engine, steps=steps, samples=samples, 
                                repeats=repeats, on_states=on_states,
                                off_states=off_states)
        res = attach_phenotype({'input_condition': sim_input, 'simul_result': res})        
        results.append(res)
    
    json.dump(results, open(file_b1, 'w'), indent=4)

# @pytest.mark.skip(reason='')
def test_summary_32results():
    with open(file_b1,'r') as fin:
        res32 = json.load(fin)

    df0 = pd.DataFrame([], columns=['input','phenotype','ratio'])
    inplabels = [
        'State_Mutagen',
        'State_GFs',
        'State_Nutrients',
        'State_TNFalpha',
        'State_Hypoxia' 
        ]
    k = 0 
    for res in res32:
        res = attach_phenotype(res)
        inpcond = res['input_condition']        
        inputstring = "".join(["1" if inpcond[x] else "0" for x in inplabels])
        attrs = res['simul_result']['attractors']
        for att in attrs: 
            phenotype = attrs[att]['phenotype']
            ratio = attrs[att]['ratio']
            df0.loc[k] = {
                'input':inputstring,
                'phenotype':phenotype,
                'ratio': ratio
                }
            k += 1
                
    df0.groupby(['input','phenotype']).sum().to_csv(file_b2)

# input 
file_a2 = join(dirname(tab_s1.__file__), 'a', 
    'a2-fumia-model-processed-weighted-sum.txt')

# output 
file_test1 = join(dirname(tab_s1.__file__), 'b', 
    'test1-simul-result-with-00000.json')
file_test2 = join(dirname(tab_s1.__file__), 'b', 
    'test2-simul-result-with-11100.json')
file_b1 = join(dirname(tab_s1.__file__), 'b', 
    'b1-simul-result-with-32-conditions.json')
file_b2 = join(dirname(tab_s1.__file__), 'b',
    'b2-simul-result-table-summary.csv')

with open(file_a2, 'r') as fobj:
    lines = fobj.readlines()
    lines2 = [] 
    for lin in lines: 
        lin = lin.strip()
        if lin[0] == '#': 
            continue 
        right = lin.split('=')[1].strip()
        if right == 'input':            
            lines2.append( lin.split('=')[0].strip() + '=' + 'False') 
        else: 
            lines2.append(lin)

    modeltext = "\n".join(lines2)

attr_cy.build(modeltext, pyx='engine.pyx', weighted_sum=True)

import pyximport; pyximport.install()

template = {
    'State_Mutagen' : False,
    'State_GFs': False,
    'State_Nutrients': False,
    'State_TNFalpha': False,
    'State_Hypoxia': False,
    'State_Gli': False,
    }

    