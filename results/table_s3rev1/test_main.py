import json, itertools, pytest, pandas as pd, pickle, sys
from os.path import exists,join,dirname
from numpy.random import random
from boolean3_addon import attr_cy, to_logic
from tqdm import tqdm 
from sbie_optdrug2.results import fumia 
from sbie_optdrug2.results.fumia import rule
from sbie_optdrug2.results import table_s3rev1

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def to_bits_in(input_data):
    
    labels = ['State_Mutagen','State_GFs','State_Nutrients','State_TNFalpha',
        'State_Hypoxia','State_Gli']        
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
        print('already exists:', chk_simul_results)
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
    for simul_input in tqdm(simul_inputs, ascii=True, 
            dynamic_ncols=True):
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


def test_postproc_1(force):

    if exists(chk_simul_results_postproc) and not force:
        print('already exists:', chk_simul_results_postproc)
        return 

    with open(chk_simul_inputs, 'rb') as f: 
        simul_inputs = pickle.load(f)

    with open(chk_simul_results, 'rb') as f: 
        simul_results = pickle.load(f)

    database = [] 

    for in0, res0 in zip(simul_inputs, simul_results):
        # print(in0['compound'])
        data0 = {} 

        data0['compound'] = in0['compound']
        data0['targetting'] = in0['targetting']
        data0['input_bits'] = res0['input_bits']
        data0['input_cond'] = in0['input_cond']
        data0['mut_config'] = res0['mut_config']
        data0['mut_bits'] = res0['mut_bits']

        attrs = res0['simul_result']['attractors']
        attr_ratio = {'A': 0.0, 'P': 0.0, 'Q': 0.0, 'U': 0.0}
        for attr_key in attrs: 
            attr_ratio[ attrs[attr_key]['phenotype'] ] += attrs[attr_key]['ratio']

        data0['A'] = attr_ratio['A']
        data0['P'] = attr_ratio['P']
        data0['Q'] = attr_ratio['Q']

        database.append(data0)

    df_simul_res = pd.DataFrame(database)

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

    assembled = df_simul_res.merge(mean_by_drug, how='left', 
        left_on='compound', right_index=True)

    assembled.to_csv(chk_simul_results_postproc, index=False)


def test_postproc_2():

    if exists(chk_simul_results_postproc_normalize): 
        return 

    d0 = pd.read_csv(chk_simul_results_postproc)

    idx_cont = d0['compound']=='control'
    idx_tr = d0['compound']!='control'    

    d0_cont = d0.loc[idx_cont]
    d0_tr = d0.loc[idx_tr]

    g_cont = d0_cont.groupby(['input_bits', 'mut_bits']).groups 
    g_tr = d0_tr.groupby(['input_bits', 'mut_bits']).groups 

    # for drug in tqdm(target_table, ascii=True): 
    for key in tqdm(g_tr, ascii=True): 
        id_cont = g_cont[key]
        id_tr = g_tr[key]
        
        Acont = d0_cont.loc[id_cont, 'A']
        Pcont = d0_cont.loc[id_cont, 'P']
        Qcont = d0_cont.loc[id_cont, 'Q']
        
        d0_tr.loc[id_tr, 'dA'] = d0_tr.loc[id_tr, 'A'] - Acont.values[0]
        d0_tr.loc[id_tr, 'dP'] = d0_tr.loc[id_tr, 'P'] - Pcont.values[0] 
        d0_tr.loc[id_tr, 'dQ'] = d0_tr.loc[id_tr, 'Q'] - Qcont.values[0]

    d0_tr.to_csv(chk_simul_results_postproc_normalize, index=False)


def test_postproc_3():
    
    if exists(output_corr_mut): 
        return 

    d0 = pd.read_csv(chk_simul_results_postproc_normalize) 


    mut_types = d0['mut_bits'].unique().tolist()
    corr_list = []
    _cols = ['ActArea','A','P','Q','dA','dP','dQ']
    for mt in mut_types: 
        dsub = d0.loc[ d0['mut_bits'] == mt ]
        meandata = dsub.groupby('compound').mean()
        rA = meandata[_cols].corr().loc['ActArea', 'A']
        rP = meandata[_cols].corr().loc['ActArea', 'P']
        rQ = meandata[_cols].corr().loc['ActArea', 'Q']
        rdA = meandata[_cols].corr().loc['ActArea', 'dA']
        rdP = meandata[_cols].corr().loc['ActArea', 'dP']
        rdQ = meandata[_cols].corr().loc['ActArea', 'dQ']
        
        corr_list.append(
            {'mut_bits': mt, 
            'A': rA, 
            'P': rP, 
            'Q': rQ, 
            'dA': rdA, 
            'dP': rdP, 
            'dQ': rdQ
            })

    pd.DataFrame(corr_list).to_csv(output_corr_mut,index=False)
    # import 

def test_postproc_4():    
    # if exists(output_corr_mut): 
    #     return 
    d0 = pd.read_csv(chk_simul_results_postproc_normalize)
    _cols = ['ActArea','A','P','Q','dA','dP','dQ']
    in_types = d0['input_bits'].unique().tolist(); corr_list = []    
    for in0 in in_types: 
        dsub = d0.loc[ d0['input_bits'] == in0 ]
        meandata = dsub.groupby('compound').mean()
        rA = meandata[_cols].corr().loc['ActArea', 'A']
        rP = meandata[_cols].corr().loc['ActArea', 'P']
        rQ = meandata[_cols].corr().loc['ActArea', 'Q']
        rdA = meandata[_cols].corr().loc['ActArea', 'dA']
        rdP = meandata[_cols].corr().loc['ActArea', 'dP']
        rdQ = meandata[_cols].corr().loc['ActArea', 'dQ']
        
        data0 = {
            'input_bits': in0, 
            'A': rA, 
            'P': rP, 
            'Q': rQ, 
            'dA': rdA, 
            'dP': rdP, 
            'dQ': rdQ
            }
        corr_list.append(data0)

    df1 = pd.DataFrame(corr_list)
    df1.to_csv(output_corr_inp,index=False)
    # import     


if not exists('fumia_engine.pyx'):
    attr_cy.build(fumia.modeltext(), 
        pyx='fumia_engine.pyx',
        weighted_sum=True
        )

import pyximport; pyximport.install()

basedir = dirname(table_s3rev1.__file__)

# dataset 
dataset_query_drugs = join(basedir,
    '../table_s2/dataset-query-drugs.csv')

dataset_ccle_dose_resp = join(basedir,
    '../table_s2/download/CCLE_NP24.2009_Drug_data_2015.02.24.csv')

dataset_ccle_crc_info = join(basedir,
    '../table_s3/dataset-ccle_crc_info.csv')

# intermediate results
chk_drug_model_target_dict = join(basedir, 
    '../table_s3/chk-drug-modeltarget.json')

chk_simul_results = join(basedir,
    'chk-simul-results.pkl')

chk_simul_inputs = join(basedir,
    'chk-simul-inputs.pkl')

chk_simul_results_postproc = join(basedir,
    'chk-simul-results-postproc.csv')

chk_simul_results_postproc_normalize = join(basedir, 
    'chk-simul-results-postproc-normal.csv')

# final results
output_corr_mut = join(basedir,
    'output-a-corr-mut.csv')

output_corr_inp = join(basedir,
    'output-b-corr-inp.csv')

max_num_drugs = 2;

samples = 10


