import json,itertools,pytest
from os.path import exists
from numpy.random import random
from boolean3_addon import attr_cy, to_logic
from pdb import set_trace
from boolean3_addon import to_logic 
from os.path import exists,dirname
import pandas as pd
from tqdm import tqdm 

from sbie_optdrug2.results import table_s1
from sbie_optdrug2.results.table_s1 import fumia_phenotype
from fumia_phenotype import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def get_input_combinations():    

    inputs = [ {
    'State_Mutagen' : State_Mutagen, 
    'State_GFs': State_GFs,
    'State_Nutrients': State_Nutrients,
    'State_TNFalpha': State_TNFalpha,
    'State_Hypoxia': State_Hypoxia,
    'State_Gli': False,
    } for State_Mutagen,State_GFs,State_Nutrients,State_TNFalpha, \
        State_Hypoxia in itertools.product([False,True], repeat=5)]

    return inputs


def get_mutations_config_all():
    list_mutations_config = []
    for State_APC, State_Ras, State_Smad, State_PTEN, State_p53 in \
        itertools.product([None, 'off', 'on'], repeat=5):

        mutations_config = {
            'State_APC': State_APC,
            'State_Ras': State_Ras,
            'State_Smad': State_Smad,
            'State_PTEN': State_PTEN,
            'State_p53': State_p53,
            }   

        list_mutations_config.append(mutations_config)

    return list_mutations_config


def get_mutations_config():
   
    # free
    mutations_config_list = [ 
        {
            'State_APC': None,
            'State_Ras': None,
            'State_Smad': None,
            'State_PTEN': None,
            'State_p53': None,
        },

        {
            'State_APC': 'off',
            'State_Ras': None,
            'State_Smad': None, 
            'State_PTEN': None, 
            'State_p53': None, 
        }, 

        {
            'State_APC': 'off',
            'State_Ras': 'on',
            'State_Smad': None, 
            'State_PTEN': None, 
            'State_p53': None, 
        }, 

        {
            'State_APC': 'off',
            'State_Ras': 'on',
            'State_Smad': 'off',
            'State_PTEN': None, 
            'State_p53': None, 
        }, 
        
        {
            'State_APC': 'off',
            'State_Ras': 'on',
            'State_Smad': 'off',
            'State_PTEN': 'off', 
            'State_p53': None, 
        },

        {
            'State_APC': 'off',
            'State_Ras': 'on',
            'State_Smad': 'off',
            'State_PTEN': 'off', 
            'State_p53': 'off',
        }
    ]
    return mutations_config_list


def test_c1(force):    
    if exists(file_c1) and force==False: 
        return

    import engine

    # mut_configs = get_mutations_config_all()
    mut_configs = get_mutations_config()
    input_combinations = get_input_combinations()

    steps = 60
    samples = 100
    repeats = 100
    ncpus = 25

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
        res = attach_phenotype({'input_condition': simul_condition, 
                                'simul_result': res})        
        results.append(res)

    json.dump(results, open(file_c1, 'w'), indent=4)


def test_c2_summary():
    with open(file_c1,'r') as fin:
        previous_results = json.load(fin)

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

    df5 = df4.sum(axis=1)

    for idx in df5.index.values: 
        df4.loc[idx] = df4.loc[idx]/df5.loc[idx]

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


# dataset 
dataset_ratio_a = join(dirname(table_s1.__file__), 'dataset-ratio-a.csv')
dataset_ratio_p = join(dirname(table_s1.__file__), 'dataset-ratio-p.csv')
dataset_ratio_q = join(dirname(table_s1.__file__), 'dataset-ratio-q.csv')

file_a2 = join(dirname(table_s1.__file__),
    'output-a1-fumia-model-processed-weighted-sum.txt')

# output 
file_c1 = join(dirname(table_s1.__file__), 'c', 
                'output-c1-results.json')

file_c2 = join(dirname(table_s1.__file__), 'c', 
                'output-c2-results-summary.csv')

file_c2_plot_a = join(dirname(table_s1.__file__), 'c', 
                    'figure-c2-simul-results-summary')

file_c2_plot_b = join(dirname(table_s1.__file__), 'c',
                    'figure-c2-simul-results-parity')

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

simul_mode = 'nohypoxia' 
# #simul_mode = 'normal'

