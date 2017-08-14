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


def test_run():
    # mut_configs = fumia.get_mutations_config()
    mut_configs = [
        {
            'State_APC': 'off',
            'State_Ras': None,
            'State_Smad': None,
            'State_PTEN': None,
            'State_p53': None,
        }
    ]

    run(mut_configs=mut_configs, targetting=[])



def run(mut_configs=[], targetting=[], 
    steps=50, samples=100, repeats=100, ncpus=50, show_progress=False):

    import fumia_engine

    input_combinations = fumia.get_input_combinations()
    input_combinations_2 = []    
    res_list = []     
    skipped_inputs = [
        '000010',
        '001010',
        '010010',
        '011010',
        '100000',
        '100010',
        '101010',
        '110010',
        '111010'
        ]
    for input0 in input_combinations:
        res = fumia.to_bits(input0)
        res_list.append(res)
        if res in skipped_inputs:
            continue 
        else: 
            input_combinations_2.append(input0)

    list_simul_condition = [] 
    for input0 in input_combinations_2:
        for mut0 in mut_configs:
            simul_condition = {'input': input0, 'mutation': mut0}
            on_states = []
            off_states = [] + targetting            
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

    for simul_condition in tqdm(list_simul_condition, ascii=True, disable=~show_progress): 
        on_states = simul_condition['on_states']
        off_states = simul_condition['off_states']
        res = attr_cy.parallel(fumia_engine, 
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
    # return results


def test_summary():
    summary() 


def summary():

    with open(chk_fumia_simul,'r') as fin:
        previous_results = json.load(fin)

    k = 0
    df0 = pd.DataFrame([], columns=['input','mutation','A','P','Q','U'])
    inplabels = ['State_Mutagen','State_GFs','State_Nutrients','State_TNFalpha','State_Hypoxia']
    mutlabels = ['State_APC','State_Ras','State_Smad','State_PTEN','State_p53']
    for res in previous_results:        
        res = rule.attach_phenotype(res)
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

    simul_mode = 'nohypoxia'
    #simul_mode = 'normal'
    # if simul_mode == 'nohypoxia':
    #     # nohypoxia가 실제상황을 더 잘 반영하는 것을 의미하는 것은 아니지만, 
    #     # 페이퍼의 결과에는 더 잘 부합한다. 
    #     df0 = df0[df0['input'] != '00001']
    #     df0 = df0[df0['input'] != '00101']
    #     df0 = df0[df0['input'] != '01001']
    #     df0 = df0[df0['input'] != '01101']
    #     df0 = df0[df0['input'] != '10000']
    #     df0 = df0[df0['input'] != '10001']
    #     df0 = df0[df0['input'] != '10101']
    #     df0 = df0[df0['input'] != '11001']
    #     df0 = df0[df0['input'] != '11101']

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
    # xxx
    df4.to_csv(chk_summary, index=False)
    
    return df4 

    
from sbie_optdrug2.results import fumia 

if not exists('fumia_engine.pyx'):
    attr_cy.build(
        fumia.modeltext(),pyx='fumia_engine.pyx', 
        weighted_sum=True
        )

import pyximport; pyximport.install()

chk_fumia_simul = 'chk-fumia-simul.json'
chk_summary = 'chk-fumia-simul-summary.csv'

