import json,itertools,pytest
from os.path import exists,join
from numpy.random import random
from boolean3_addon import attr_cy, to_logic
from pdb import set_trace
from boolean3_addon import to_logic 
from os.path import exists,dirname
import pandas as pd
from tqdm import tqdm 

from sbie_optdrug2.results import fumia 
from sbie_optdrug2.results.fumia import rule

# from sbie_optdrug.result_season2 import tab_s1
# from sbie_optdrug.result_season2.tab_s1 import fumia_phenotype
# from fumia_phenotype import *

# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt


def test_c1():    

    import engine

    mut_configs = fumia.get_mutations_config()
    input_combinations = fumia.get_input_combinations()

    steps = 60
    samples = 100
    repeats = 100
    ncpus = 60

    list_simul_condition = [] 
    for input0 in input_combinations:
        for mut0 in mut_configs:
            simul_condition = {'input': input0, 'mutation': mut0}
            on_states = [] 
            off_states = [] 

            for lbl in input0: 
                if input0[lbl] == True: 
                    on_states.append(lbl)
                else:
                    off_states.append(lbl)

            for lbl in mut0: 
                if mut0[lbl] == 'on': 
                    on_states.append(lbl)
                elif mut0[lbl] == 'off':
                    off_states.append(lbl)
                else: 
                    pass

            list_simul_condition.append({
                'simul_condition': simul_condition,
                'on_states': on_states, 
                'off_states': off_states, 
                })

    results = []

    for simul_condition in tqdm(list_simul_condition, ascii=True): 
        on_states = simul_condition['on_states']
        off_states = simul_condition['off_states']
        res = attr_cy.parallel(engine, 
                                steps=steps,
                                samples=samples,
                                ncpus=ncpus,
                                repeats=repeats, 
                                on_states=on_states,
                                off_states=off_states
                                )
        res = rule.attach_phenotype({'input_condition': simul_condition, 
                                'simul_result': res})        
        results.append(res)

    json.dump(results, open(chk_fumia_simul, 'w'), indent=4)


@pytest.mark.skip(reason='')
def test_c2_summary():

    with open(chk_fumia_simul,'r') as fin:
        previous_results = json.load(fin)

    simul_mode = 'nohypoxia'
    #simul_mode = 'normal'

    df0 = pd.DataFrame([], columns=['input','mutation','A','P','Q','U'])

    k = 0 

    for res in previous_results:        
        res = attach_phenotype(res)
        inpcond = res['input_condition']['simul_condition']['input']
        inputstring = "".join(["1" if inpcond[x] else "0" for x in inplabels])
        mutcond = res['input_condition']['simul_condition']['mutation']
        mutstring = "".join(["0" if mutcond[x] == None else "1" for x in mutlabels])
        attrs = res['simul_result']['attractors']
        df0.loc[k, 'input'] = inputstring
        df0.loc[k, 'mutation'] = mutstring    
        df0.loc[k, 'A'] = 0.0
        df0.loc[k, 'P'] = 0.0
        df0.loc[k, 'Q'] = 0.0
        df0.loc[k, 'U'] = 0.0
        for att in attrs: 
            df0.loc[k, attrs[att]['phenotype']] += attrs[att]['ratio']

        k += 1

    if simul_mode == 'nohypoxia':
        df0 = df0[df0['input'] != '00001']
        df0 = df0[df0['input'] != '00101']
        df0 = df0[df0['input'] != '01001']
        df0 = df0[df0['input'] != '01101']
        df0 = df0[df0['input'] != '10000']
        df0 = df0[df0['input'] != '10001']
        df0 = df0[df0['input'] != '10101']
        df0 = df0[df0['input'] != '11001']
        df0 = df0[df0['input'] != '11101']

    df2 = df0.groupby(['mutation','input']).sum()
    df3 = df2.reset_index(level=df2.index.names)
    
    df4 = df3.groupby('mutation').sum()
    # (Pdb++) df4
    #                  A          P          Q        U
    # mutation
    # 00000     13.85450   0.501875   1.643625  0.00000
    # 10000     27.70900   2.003750   2.287250  0.00000
    # 11000     27.57100   3.008000   1.421000  0.00000
    # 11100     25.32875   5.018500   1.636500  0.01625
    # 11110     19.19550   6.455000   6.349500  0.00000
    # 11111     11.53900  10.000000  10.461000  0.00000

    df5 = df4.sum(axis=1)
    # (Pdb++) df5
    # mutation
    # 00000    16.0
    # 10000    32.0
    # 11000    32.0
    # 11100    32.0
    # 11110    32.0
    # 11111    32.0
    # dtype: float64

    for idx in df5.index.values: 
        df4.loc[idx] = df4.loc[idx]/df5.loc[idx]

    # ipdb> df4
    #                  A         P         Q         U
    # mutation
    # 00000     0.865906  0.031367  0.102727  0.000000
    # 10000     0.865906  0.062617  0.071477  0.000000
    # 11000     0.861594  0.094000  0.044406  0.000000
    # 11100     0.791523  0.156828  0.051141  0.000508
    # 11110     0.599859  0.201719  0.198422  0.000000
    # 11111     0.360594  0.312500  0.326906  0.000000

    plt.figure(); 
    df4.plot(marker='o');     
    plt.savefig(file_c2_plot_a + '-' + simul_mode)

    fumia_a = pd.read_csv('a/a4-ratio-a.csv',index_col='x')
    fumia_p = pd.read_csv('a/a4-ratio-p.csv',index_col='x')
    fumia_q = pd.read_csv('a/a4-ratio-q.csv',index_col='x')

    plt.figure(); 
    line1, = plt.plot(df4['A'], fumia_a['Curve1'],'o', label='A');
    line2, = plt.plot(df4['P'], fumia_p['Curve1'], 'o', label='P');
    line3, = plt.plot(df4['Q'], fumia_q['Curve1'],'o', label='Q');

    plt.plot([0,1],[0,1], color='0.75', alpha=0.4);
    plt.legend(handles=[line1,line2,line3])    
    plt.xlabel('Reproduced data')
    plt.ylabel('Original data')    
    plt.savefig(file_c2_plot_b + '-' + simul_mode)

# # input 
# file_a2 = join(dirname(tab_s1.__file__), 'a', 
#                 'a2-fumia-model-processed-weighted-sum.txt')

# # output 
# file_c1 = join(dirname(tab_s1.__file__), 'c', 
#                 'c1-simul-results.json')

# file_c2 = join(dirname(tab_s1.__file__), 'c', 
#                 'c2-simul-results-summary.csv')

# file_c2_plot_a = join(dirname(tab_s1.__file__), 'c', 
#                     'c2-simul-results-summary')

# file_c2_plot_b = join(dirname(tab_s1.__file__), 'c',
#                     'c2-simul-results-parity')


# with open(file_a2, 'r') as fobj:
#     lines = fobj.readlines()
#     lines2 = [] 
#     for lin in lines: 
#         lin = lin.strip()
#         if lin[0] == '#': 
#             continue 
#         right = lin.split('=')[1].strip()
#         if right == 'input':            
#             lines2.append( lin.split('=')[0].strip() + '=' + 'False') 
#         else: 
#             lines2.append(lin)

#     modeltext = "\n".join(lines2)

from sbie_optdrug2.results import fumia 

attr_cy.build(
    fumia.modeltext(),
    pyx='engine.pyx',
    weighted_sum=False
    )

import pyximport; pyximport.install()

inplabels = [
    'State_Mutagen',
    'State_GFs',
    'State_Nutrients',
    'State_TNFalpha',
    'State_Hypoxia' 
    ]

mutlabels = [
    'State_APC',
    'State_Ras',
    'State_Smad',
    'State_PTEN',
    'State_p53' 
    ]

chk_fumia_simul = 'chk-fumia-simul.json'