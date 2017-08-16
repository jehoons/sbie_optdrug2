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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle

def to_bits_in(input_data):
    labels = ['State_Mutagen','State_GFs','State_Nutrients','State_TNFalpha','State_Hypoxia','State_Gli']
    res = ['1' if input_data[state] else '0' for state in labels]
    return 'b'+"".join(res)

def to_bits_mut(input_data):    
    labels = ['State_APC','State_Ras','State_Smad','State_PTEN','State_p53']
    res = ['1' if input_data[state] else '0' for state in labels]
    return 'b'+"".join(res)    

def run(input_cond=None, mut_config=None, targetting=[], 
    steps=50, samples=100, repeats=100, ncpus=64):
    import fumia_engine
    on_states = []
    off_states = [] + targetting
    for lbl in input_cond: 
        if input_cond[lbl] == True: 
            on_states.append(lbl)
        else:
            off_states.append(lbl)

    for lbl in mut_config: 
        if mut_config[lbl] == 'on': 
            on_states.append(lbl)
        elif mut_config[lbl] == 'off':
            off_states.append(lbl)

    res = attr_cy.parallel(fumia_engine, 
                            steps=steps,
                            samples=samples,
                            ncpus=ncpus,
                            repeats=repeats, 
                            on_states=on_states,
                            off_states=off_states
                            )
    res = rule.attach_phenotype(
            {
                'input_condition': input_cond, 
                'input_bits': to_bits_in(input_cond),
                'mut_config': mut_config, 
                'mut_bits': to_bits_mut(mut_config),
                'targetting': targetting, 
                'simul_result': res
            }
        )

    return res


def test_run_simul(force):

    if exists(chk_simul_results) and force == False:
        return 

    with open(chk_drug_model_target_dict, 'r') as f: 
        target_table = json.load(f)

    mut_configs = fumia.get_mutations_config()
    input_configs = fumia.get_input_combinations()

    target_table['control'] = [] 
    pert_res = {}

    simul_inputs = [] 
    for mut0 in mut_configs:
        for input0 in input_configs:
            # for drug in tqdm(target_table, ascii=True): 
            for drug in target_table:
                targets = target_table[drug]
                if len(targets)>= max_num_drugs: 
                    targets = targets[0:max_num_drugs]

                simul_input_0 = {
                    'input_cond': input0, 
                    'mut_config': mut0,
                    'targetting': targets, 
                    'compound': drug
                    }

                simul_inputs.append(simul_input_0)
    
    simul_results = [] 
    for simul_input in tqdm(simul_inputs, ascii=True):
        res = run(input_cond=simul_input['input_cond'],
            mut_config=simul_input['mut_config'],
            targetting=simul_input['targetting'],
            samples=samples
            )

        simul_results.append(res)

    with open(chk_simul_results, 'wb') as f: 
        pickle.dump(simul_results, f)

    with open(chk_simul_inputs, 'wb') as f: 
        pickle.dump(simul_inputs, f)        


def test_postproc():

    with open(chk_simul_inputs, 'rb') as f: 
        simul_inputs = pickle.load(f)

    with open(chk_simul_results, 'rb') as f: 
        simul_results = pickle.load(f)

    database = [] 

    for in0, res0 in zip(simul_inputs, simul_results):

        data0 = {} 

        data0['compound'] = in0['compound']
        data0['targetting'] = in0['targetting']
        data0['input_bits'] = res0['input_bits']
        data0['input_cond'] = in0['input_cond']
        data0['mut_config'] = res0['mut_config']
        data0['mut_bits'] = res0['mut_bits']

        attrs = res0['simul_result']['attractors']
        attr_ratio = {'A': 0, 'P': 0, 'Q': 0}
        for attr_key in attrs: 
            attr_ratio[ attrs[attr_key]['phenotype'] ] += attrs[attr_key]['ratio']

        data0['A'] = attr_ratio['A']
        data0['P'] = attr_ratio['P']
        data0['Q'] = attr_ratio['Q']

        database.append(data0)

    df_simul_res = pd.DataFrame(database)

    # pert_res = pd.read_json(chk_pert_res).transpose()
    query_drugs = pd.read_csv(dataset_query_drugs)
    ccle_dose_resp = pd.read_csv(dataset_ccle_dose_resp)
    ccle_crc = pd.read_csv(dataset_ccle_crc_info)
    cells_df = ccle_crc[['Cell line primary name']].drop_duplicates()
    
    ccle_dose_resp_sel = ccle_dose_resp.merge(
        cells_df,
        how='inner',
        left_on='Primary Cell Line Name',
        right_on='Cell line primary name'
        )
    
    ccle_dose_resp_sel = ccle_dose_resp
    ccle_drug_chemblid = query_drugs[['compound_name', 
            'compound_chembl_id']].drop_duplicates()
    
    # # with compound_chembl_id 
    ccle_dose_resp_sel_2 = ccle_dose_resp_sel.merge(
        ccle_drug_chemblid,
        how='inner',
        left_on='Compound',
        right_on='compound_name'
        )

    mean_by_drug = ccle_dose_resp_sel_2.groupby('compound_chembl_id').mean()

    assembled = df_simul_res.merge(
        mean_by_drug, 
        how='inner',
        left_on='compound',
        right_index=True
        )

    assembled.to_csv(chk_simul_results_postproc)


if not exists('fumia_engine.pyx'):
    attr_cy.build(
        fumia.modeltext(),pyx='fumia_engine.pyx', 
        weighted_sum=True
        )

import pyximport; pyximport.install()

dataset_ccle_crc_info = 'dataset-ccle_crc_info.csv'
dataset_query_drugs = '../table_s2/dataset-query-drugs.csv'
dataset_ccle_dose_resp = '../table_s2/download/CCLE_NP24.2009_Drug_data_2015.02.24.csv'

chk_drug_model_target_dict = '../table_s3/chk-drug-modeltarget.json'
chk_simul_results = 'chk-simul-results.pkl'
chk_simul_inputs = 'chk-simul-inputs.pkl'
chk_simul_results_postproc = 'chk-simul-results-postproc.csv'

max_num_drugs = 2;
samples = 1000
