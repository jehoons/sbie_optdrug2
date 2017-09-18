import json, itertools, pytest, pandas as pd, pickle, sys
from os.path import exists,join,dirname
from numpy.random import random
from boolean3_addon import attr_cy, to_logic
from tqdm import tqdm 
from sbie_optdrug2.results import table_s3rev2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import model_fumia2013
from model_fumia2013 import rule 
from ipdb import set_trace
import dataset_msigdb
from pandas.io.json import json_normalize
import seaborn as sns, numpy as np

import time
import numpy as np
from sklearn import datasets, linear_model    
from sklearn.metrics import mean_squared_error, r2_score

basedir = dirname(table_s3rev2.__file__)

# Dataset 
dataset_query_drugs = join(basedir,
    '../table_s2/dataset-query-drugs.csv')

dataset_ccle_dose_resp = join(basedir,
    '../table_s2/download/CCLE_NP24.2009_Drug_data_2015.02.24.csv')

# CRC subset of CCLE 
dataset_ccle_crc_info = join(basedir,
    '../table_s3/dataset-ccle_crc_info.csv')

chk_drug_model_target_dict = join(basedir, 
    '../table_s3/chk-drug-modeltarget.json')

msigdb = dataset_msigdb.load_msigdb()

if not exists('fumia_engine.pyx'):
    attr_cy.build(model_fumia2013.modeltext(), 
        pyx='fumia_engine.pyx',
        weighted_sum=True
        )

import pyximport; pyximport.install()

max_num_drugs = 2;

samples = 10


with open('dataset-goterms-in-apq.json', 'r') as f: 
    goterms_in_apq = json.load(f)


def genesymdict():

    fum_gmap = model_fumia2013.genesym_map()
    fum_gmap.set_index('From', inplace=True)
    
    fum_nodeinfo = model_fumia2013.nodeinfo()
    fum_nodeinfo['uniprot_id']=fum_nodeinfo['uniprot_id'].astype(str)

    node_uniprot_map = {} 
    for i in fum_nodeinfo.index:
        varname = fum_nodeinfo.loc[i, 'variable']
        uniprot_id = fum_nodeinfo.loc[i, 'uniprot_id'].split(';')
        if uniprot_id == ['nan']: 
            uniprot_id = None
        node_uniprot_map[varname] = uniprot_id

    node_to_genesym = {} 

    for node in node_uniprot_map:         
        if node_uniprot_map[node] == None: 
            node_to_genesym[node] = [] 
            continue 

        gslist = [] 
        for uni in node_uniprot_map[node]:
            if uni in fum_gmap.index:
                gslist.append(fum_gmap.loc[uni, 'To'])

        node_to_genesym[node] = gslist

    d = dict((k,v) for k, v in node_to_genesym.items() if v!=[])
    return d


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


def mygsea(on_statestates_genesyms):
    
    ''' 켜져 있는 상태의 유전자는 어느 GO에 속하는가? '''
    
    if exists('chk-fumia-msigdb.csv'): 
        return 
    
    import dataset_msigdb
    
    msigdb = dataset_msigdb.load_msigdb()
    fumia_msigdb = {} 
    fumia_msigdb_list = [] 
    for k in msigdb.keys(): 
        mutual_genes = set(msigdb[k]['geneset']).intersection(set(on_states_genesyms))
        data = {
            'signature': k, 
            'mutual_genes': [x for x in mutual_genes],
            'num_interection': len(mutual_genes),
            'ratio_interection': len(mutual_genes) / len(msigdb[k]['geneset']),
            # 'geneset': msigdb[k]['geneset'],
            # 'on_states_genesyms': on_states_genesyms, 
            # 'on_states': on_states
            }
        fumia_msigdb[k] = data 
        fumia_msigdb_list.append(data) 

    with open('chk-fumia-msigdb.json', 'w') as f: 
        json.dump(fumia_msigdb,f,indent=4)

    df = pd.DataFrame(fumia_msigdb_list)
    df.to_csv('chk-fumia-msigdb.csv') 


def mygsea_rev2(on_states_genesyms):

    ''' 켜져 있는 상태의 유전자는 어느 클래스의 GO에 속하는가? '''
    
    apq_sel = goterms_in_apq['selected']

    data_list = []     
    for cl in apq_sel:
        for goterm in apq_sel[cl]:
            mutual_genes = set(msigdb[goterm]['geneset']).intersection(set(on_states_genesyms))
            data = {
                'class': cl, 
                'signature': goterm, 
                'mutual_genes': [x for x in mutual_genes],
                'num_interection': len(mutual_genes),
                'geneset_size': len(msigdb[goterm]['geneset']), 
                'ratio_interection': len(mutual_genes) / len(msigdb[goterm]['geneset']),
                }
            data_list.append(data)

    return data_list 


def get_on_states_genes(labels, statekey):
    
    ''' 여기서 msigdb는 uniport_id로 되어있다. 그러므로 gene symbol으로 바꾸어줄 필요가 있다.
    1. uniprot to gene symbol '''

    id_list = [] 
    
    import numpy as np 
    
    for id in model_fumia2013.nodeinfo()['uniprot_id'].unique():
        if type(id) == str: 
            id_list += id.split(';')
    
    pd.DataFrame(id_list).to_csv('chk-fumia-uniprot.csv', index=False)
            
    fum_gmap = model_fumia2013.genesym_map()

    on_states = [ l if s=='1' else None for s,l in zip(statekey, labels)]
    on_states = set(on_states)
    on_states.remove(None)

    fum_nodeinfo = model_fumia2013.nodeinfo()
    fum_nodeinfo['uniprot_id']=fum_nodeinfo['uniprot_id'].astype(str)

    node_uniprot_map = {} 
    for i in fum_nodeinfo.index:
        varname = fum_nodeinfo.loc[i, 'variable']
        uniprot_id = fum_nodeinfo.loc[i, 'uniprot_id'].split(';')
        if uniprot_id == ['nan']: 
            uniprot_id = None
        node_uniprot_map[varname] = uniprot_id
    
    fumdict = fum_gmap.to_dict(orient='list')
    fumdict2 = {} 
    for to0, frm in zip(fumdict['To'], fumdict['From']):
        fumdict2[frm] = to0

    grp = fum_nodeinfo.groupby('variable').groups

    on_states_uniprots = [] 

    for state in on_states: 
        ids = grp[state]
        on_states_uniprots += fum_nodeinfo.loc[ids, 'uniprot_id'].tolist()

    on_states_uniprots2 = [] 
    for uni in on_states_uniprots: 
        uni = uni.split(';')
        on_states_uniprots2 += uni 

    on_states_uniprots2 = set(on_states_uniprots2)
    on_states_uniprots2.remove('nan')

    on_states_genesyms = [] 
    for x in on_states_uniprots2: 
        if x in fumdict2:
            gs = fumdict2[x]
            on_states_genesyms.append(gs) 
    
    return on_states_genesyms  


# test: test_run_simul
# output: 
chk_simul_results = join(basedir, 'chk-simul-results.pkl')
chk_simul_inputs = join(basedir, 'chk-simul-inputs.pkl')
def test_run_simul(force): 

    if exists(chk_simul_inputs) and force == False:
        return                
    with open(chk_drug_model_target_dict, 'r') as f: 
        target_table = json.load(f)

    mut_configs = model_fumia2013.get_mutations_config()
    input_configs = model_fumia2013.get_input_combinations()

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


chk_simul_results_postproc_gsea = join(basedir,
    'chk-simul-results-postproc-gsea.csv')
def test_postproc_1(force):

    if exists(chk_simul_results_postproc_gsea) and not force:
        return

    with open(chk_simul_inputs, 'rb') as f: 
        simul_inputs = pickle.load(f)

    with open(chk_simul_results, 'rb') as f: 
        simul_results = pickle.load(f)

    idx = 1; res2 = []
    for sim_result_1 in tqdm(simul_results):
        sim_result_1 = sim_result_1['simul_result']
        labels = sim_result_1['labels']
        for attr_key in sim_result_1['attractors']: 
            state_key_values = sim_result_1['attractors'][attr_key]['value']
            ratio = sim_result_1['attractors'][attr_key]['ratio']
            phenotype = sim_result_1['attractors'][attr_key]['phenotype']
            attr_type = sim_result_1['attractors'][attr_key]['type']

            if type(state_key_values) != list: 
                state_key_values = [state_key_values] 
            for state_key in state_key_values: 
                bits = sim_result_1['state_key'][state_key]
                labels = sim_result_1['labels']            
                on_states_genesyms = get_on_states_genes(   
                    labels, bits)
                res = mygsea_rev2(on_states_genesyms)                
                for res0 in res: 
                    res0['simul_id'] = idx
                    res0['attr_key'] = attr_key
                    res0['state_key'] = state_key
                    res0['attr_ratio'] = ratio 
                    res0['attr_phenotype'] = phenotype
                    res0['attr_type'] = attr_type
                    # set_trace()
                    res2.append(res0)
        
        idx += 1

    pd.DataFrame(res2).to_csv(chk_simul_results_postproc_gsea)


# test: test_postproc_make_df
# outputs: 
chk_simul_res_rec = join(basedir,
        'chk-simul-res-rec.csv')
chk_simul_res_state_rec = join(basedir,
        'chk-simul-res-state-rec.csv')
chk_simul_res_stategenesymb_rec = join(basedir, 
        'chk-simul-res-stategenesymb-rec.csv')
def test_postproc_make_df(force):
    ''' 여기서는 시뮬레이션 입력데이터와 결과데이터를 통합하여 CSV 형식의 
    자료를 만든다. '''
    if exists('chk-simul-res-stategenesymb-rec.csv') and not force:
        return

    with open(chk_simul_inputs, 'rb') as f: 
        simul_inputs = pickle.load(f)

    with open(chk_simul_results, 'rb') as f: 
        simul_results = pickle.load(f)

    genedict = genesymdict()

    simul_results_list = [] 
    simul_results_state_list = [] 
    simul_results_stategenesymb_list = []
    idx = 1
    for inp, res in zip(simul_inputs, simul_results):
        att = res['simul_result']['attractors']
        for akey in att: 
            value = att[akey]['value']
            fumia_phenotype = att[akey]['phenotype']
            if type(value) != list: 
                value = [value]
            attr_type = att[akey]['type']
            ratio = att[akey]['ratio']

            for state_key in value:                
                flat_data = {
                    'simul_id': idx, 
                    'input_cond': inp['input_cond'], 
                    'mut_config': inp['mut_config'],
                    'compound': inp['compound'],
                    'targetting': ";".join(inp['targetting']),
                    'attr_key': akey,
                    'state_key': state_key,
                    'attr_ratio': ratio,
                    'attr_ratio_2': ratio/len(value),
                    'attr_type': attr_type,
                    'fumia_phenotype': fumia_phenotype, 
                    }

                simul_results_list.append(flat_data)
                statebits = res['simul_result']['state_key'][state_key]
                lbl = res['simul_result']['labels']                
                flat_data2 = {} 
                for l,s in zip(lbl,statebits):
                    flat_data2[l] = 1 if s=='1' else 0 
                
                simul_results_state_list.append(flat_data2)

                flat_data3 = {}
                for k in flat_data2:
                    if k in genedict: 
                        genesyms = genedict[k]
                        state_value = flat_data2[k]
                        for g in genesyms: 
                            flat_data3[g] = state_value

                simul_results_stategenesymb_list.append(flat_data3)
        
        idx += 1 

    json_normalize(simul_results_list).to_csv(
        chk_simul_res_rec, index=False)
    
    pd.DataFrame(simul_results_state_list).to_csv(
        chk_simul_res_state_rec, index=False)
    
    pd.DataFrame(simul_results_stategenesymb_list).to_csv(
        chk_simul_res_stategenesymb_rec, index=False)


# test: test_postproc_combine_simul_and_gsea
# outputs: 
output_combined_with_gsea = join(basedir, 
    'output-combined-with-gsea.csv')
def test_postproc_combine_simul_and_gsea(force):

    if exists(output_combined_with_gsea) and not force: 
        return 

    with open(chk_simul_inputs, 'rb') as f:
        simul_inputs = pickle.load(f)

    df1 = pd.read_csv('chk-simul-res-rec.csv')
    df2 = pd.read_csv(chk_simul_results_postproc_gsea)

    rois = [
        'compound',
        'input_cond.State_GFs',
        'input_cond.State_Gli',
        'input_cond.State_Hypoxia',
        'input_cond.State_Mutagen',
        'input_cond.State_Nutrients',
        'input_cond.State_TNFalpha',
        'mut_config.State_APC',
        'mut_config.State_PTEN',
        'mut_config.State_Ras',
        'mut_config.State_Smad',
        'mut_config.State_p53',
        'simul_id',
        'targetting'
        ] 

    df_simul_cond = df1.groupby('simul_id')[rois].first().copy()
    df3 = df2.merge(df_simul_cond, how='left', left_on='simul_id', right_on='simul_id')
    df3.to_csv(output_combined_with_gsea, index=False)


# test: test_plot 
# outputs: 
output_simul_result_mean_score = join(basedir,
        'output-simul-result-mean-score.csv')
fig_predicted_score_hist = 'fig-predicted-score-hist.png'
fig_predicted_score_sorted = 'fig-predicted-score-sorted.png'
def test_plot(force):

    if exists(output_simul_result_mean_score) and not force: 
        return 

    df0 = pd.read_csv(output_combined_with_gsea)
    df1 = df0[['ratio_interection','class']].copy()

    grps = df1.groupby('class').groups

    v_mean = df1.groupby('class').mean()
    v_std = df1.groupby('class').std()

    for key in grps: 
        ids = grps[key]
        v_mean0 = v_mean.loc[key].tolist()[0]
        v_std0 = v_std.loc[key].tolist()[0]
        df1.loc[ids, 'ratio_interection'] = (
            df1.loc[ids, 'ratio_interection']-v_mean0)/v_std0

    plt.figure()
    
    g = sns.FacetGrid(df1,  col="class", margin_titles=True)
    bins = np.arange(df1['ratio_interection'].min(),
            df1['ratio_interection'].max(),0.05)    
    g.map(plt.hist, 'ratio_interection')
    
    plt.savefig(fig_predicted_score_hist)

    df2 = df0[['compound', 'ratio_interection','class']].copy()

    df_arrest = df2[df2['class']=='ARREST'].groupby('compound').mean()
    df_apoptosis = df2[df2['class']=='APOPTOSIS'].groupby('compound').mean()
    df_proliferation = df2[df2['class']=='PROLIFERATION'].groupby('compound').mean()

    df_mean1 = pd.DataFrame([], columns=['compound','ratio_interection','class'])
    df_mean2 = pd.DataFrame([], columns=['compound','ratio_interection','class'])
    df_mean3 = pd.DataFrame([], columns=['compound','ratio_interection','class'])    

    df_mean1['compound'] = df_arrest.index.tolist()
    df_mean1['ratio_interection'] = df_arrest['ratio_interection'].tolist() 
    df_mean1['class'] = ['ARREST'] * df_arrest.shape[0]

    df_mean2['compound'] = df_apoptosis.index.tolist()
    df_mean2['ratio_interection'] = df_apoptosis['ratio_interection'].tolist() 
    df_mean2['class'] = ['APOPTOSIS'] * df_apoptosis.shape[0]

    df_mean3['compound'] = df_proliferation.index.tolist()
    df_mean3['ratio_interection'] = df_proliferation['ratio_interection'].tolist() 
    df_mean3['class'] = ['PROLIFERATION'] * df_proliferation.shape[0]

    df_final = pd.concat([df_mean1, df_mean2, df_mean3], ignore_index=True)

    df_final_sorted = df_final.sort_values(by='ratio_interection')
    vv = df_final_sorted['ratio_interection']
    std_v = df_final_sorted['ratio_interection'].std()
    mean_v = df_final_sorted['ratio_interection'].mean()
    df_final_sorted['ratio_interection'] = (vv - mean_v) / std_v
        
    grps = df_final_sorted.groupby('class').groups
    for key in grps: 
        ids = grps[key]
        mu0 = df_final_sorted.loc[ids, 'ratio_interection'].mean()
        std0 = df_final_sorted.loc[ids, 'ratio_interection'].std()
        df_final_sorted.loc[ids, 'ratio_interection'] = (
            df_final_sorted.loc[ids, 'ratio_interection']-mu0)/std0

    plt.figure()

    g = sns.FacetGrid(df_final_sorted,  col="class",    
        margin_titles=True, col_order=['PROLIFERATION', 'ARREST', 'APOPTOSIS'])
    
    g.map(sns.barplot, 'compound', 'ratio_interection')

    plt.savefig(fig_predicted_score_sorted)

    df_final_sorted.to_csv(output_simul_result_mean_score)


# test: test_expr_vs_predict
# outputs: 
chk_merged_reduced_mean = 'chk-merged-reduced-mean.csv'
fig_pairplot = 'fig-pairplot.png'
fig_pairplot_mean = 'fig-pairplot-mean.png'
def test_expr_vs_predict(force):

    if exists(chk_merged_reduced_mean) and not force: 
        return 

    dfcell = pd.read_csv(
            '../table_s3/dataset-ccle_crc_info.csv')
    dfdrugs = pd.read_csv(
            '../table_s2/dataset-query-drugs.csv')
    dftreat = pd.read_csv(
            '../table_s2/download/CCLE_NP24.2009_Drug_data_2015.02.24.csv')

    dfdrugs = dfdrugs[['compound_name','compound_chembl_id']
        ].drop_duplicates()

    dfexpr = dftreat.merge(
            dfcell,
            how='left',
            left_on='Primary Cell Line Name',
            right_on='Cell line primary name'
            )

    dfexpr = dfexpr.dropna(axis=0, subset=['Histology'])

    dfexpr = dfexpr.merge(
            dfdrugs, 
            how='left', 
            left_on='Compound', 
            right_on='compound_name'
            )
    
    dfsimulres = pd.read_csv(
            'output-simul-result-mean-score.csv'
            )[['compound','ratio_interection','class']]

    dfsimulres_grps = dfsimulres.groupby('class').groups    
    dfsimulres_class_list = [] 
    for k in dfsimulres_grps:         
        _df = dfsimulres.loc[dfsimulres_grps[k],
            ['compound','ratio_interection']]
        _df = _df[ _df['compound'] != 'control' ]
        _df.rename(columns={'ratio_interection': k}, inplace=True)
        _df.set_index('compound', inplace=True)
        dfsimulres_class_list.append(_df)

    dfsimulres_transp = pd.concat(dfsimulres_class_list, axis=1)

    df_merged = dfexpr.merge(dfsimulres_transp, how='left', 
            left_on='compound_chembl_id', right_index=True)
    
    df_merged_reduced = df_merged[['Primary Cell Line Name','APOPTOSIS', 
        'ARREST','PROLIFERATION','ActArea']]
    df_merged_reduced = df_merged_reduced.dropna()
    
    plt.figure() 

    g = sns.pairplot(df_merged_reduced, hue='Primary Cell Line Name', kind='reg')

    plt.savefig(fig_pairplot)

    df_merged_reduced = df_merged[['compound_chembl_id','APOPTOSIS', 
        'ARREST','PROLIFERATION','ActArea']]
    df_merged_reduced_mean = df_merged_reduced.groupby('compound_chembl_id').mean()
    df_merged_reduced_mean = df_merged_reduced_mean.dropna()

    plt.figure() 
    
    g = sns.pairplot(df_merged_reduced_mean, kind='reg')
    
    plt.savefig(fig_pairplot_mean)
    
    df_merged_reduced_mean.to_csv(chk_merged_reduced_mean)


# test: test_regression
# outputs: 
fig_reg_parity = 'fig-reg-parity.png'
def test_regression(force):

    # if exists(fig_reg_parity) and not force:
    #     return 

    fumiacrc = pd.read_csv(chk_merged_reduced_mean) 
    X = fumiacrc[['APOPTOSIS','ARREST','PROLIFERATION']].values
    # X = fumiacrc[['ARREST']].values
    y = fumiacrc['ActArea'].values

    regr = linear_model.LinearRegression()
    regr.fit(X, y)
    y_pred = regr.predict(X)

    # The coefficients
    print('Coefficients: \n', regr.coef_)
    # The mean squared error
    print("Mean squared error: %.2f" % mean_squared_error(y, y_pred))
    # Explained variance score: 1 is perfect prediction
    print('Variance score: %.2f' % r2_score(y, y_pred))

    plt.figure(figsize=(5,5))

    plt.scatter(y, y_pred, color='black') 
    y_max = max(y.tolist()+y_pred.tolist())
    y_min = min(y.tolist()+y_pred.tolist())    
    plt.plot([y_min, y_max],[y_min, y_max])
    plt.xlabel('ActArea(Experiment)'); 
    plt.ylabel('ActArea(Prediction)')
    
    plt.savefig(fig_reg_parity) 
