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
from sbie_optdrug2.results.fumia import simulator

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def test_simulate_with_ccle_drugs(with_small, force):

    if exists(chk_pert_res) and not force:
        print('file exists: ', chk_pert_res)
        return 

    dataset_fumia_node_info = fumia.dataset_fumia_node_info 
    dataset_query_drugs = '../table_s2/dataset-query-drugs.csv'

    df_model_info = pd.read_csv(dataset_fumia_node_info) 
    df_query = pd.read_csv(dataset_query_drugs) 

    with open('../table_s2/output-a-alternative-targets.json') as f:
        jsondata = json.load(f)

    query_drugs = jsondata['similarity>=0.35'].keys() 
    fumia_nodes = df_model_info['chembl_id'].dropna().unique().tolist()
    found_dict = {}
    targetting_table = pd.DataFrame([], columns=['compound', 'target', 'model_variable', 'score'])
    i = 0 
    for q0 in query_drugs:
        conn_data = jsondata['similarity>=0.35'][q0]['connectivity_w']
        targets = [t for t in conn_data.keys()]
        found_dict[q0] = 0
        for t0 in targets:
            score = conn_data = jsondata['similarity>=0.35'][q0]['connectivity_w'][t0]
            if t0 in fumia_nodes:
                found_dict[q0] += 1
                model_nodes = df_model_info.loc[
                    df_model_info['chembl_id'] == t0, 'variable'].tolist()
                for node in model_nodes: 
                    targetting_table.loc[i, 'compound'] = q0
                    targetting_table.loc[i, 'target'] = t0
                    targetting_table.loc[i, 'model_variable'] = node
                    targetting_table.loc[i, 'score'] = float(score)
                    i+=1 

    targetting_table.to_csv(chk_target_table)
    targetting_table['score'] = targetting_table['score'].astype(float)
        
    def returnThresh(x):
        sorted0 = targetting_table_updated.loc[x.index].sort_values(
            by='score', ascending=False) 
        return ";".join( sorted0['model_variable'].tolist() )

    def returnBest(x):
        score_max = targetting_table.loc[x.index, 'score'].max()
        idx_max = targetting_table['score'] == score_max
        variable_best = targetting_table.loc[idx_max, 'model_variable'].tolist()[0]
        return variable_best

    targetting_table_updated = targetting_table.groupby(['compound','model_variable']).mean().reset_index()
    targetting_table_updated2 = targetting_table_updated.groupby('compound')['model_variable'].apply(lambda x: returnThresh(x))
    
    targetting_table_updated2_json = {} 

    for com in targetting_table_updated2.index: 
        var0 = targetting_table_updated2[com]
        targetting_table_updated2_json[com] = var0.split(';')


    with open(chk_drug_model_target_dict, 'w') as f: 
        json.dump(targetting_table_updated2_json, f, indent=4)

    targetting_table_updated2_json['control'] = [] 
    pert_res = {}
    for drug in tqdm(targetting_table_updated2_json, ascii=True): 
        targets = targetting_table_updated2_json[drug]

        # max_num_drugs = 2;
        if len(targets)>= max_num_drugs: targets = targets[0:max_num_drugs]
        
        simulator.run(mut_configs=mut_configs, targetting=targets, samples=10)
        res = simulator.summary()        
        pert_res[drug] = res.to_dict(orient='records')[0]    
        # break 

    # xxx
    with open(chk_pert_res, 'w') as f: 
        json.dump(pert_res, f, indent=4)


def test_postproc_simulation_results(with_small, force):
    if exists(chk_fumia_with_ccle_data_cell_mean) and force==False:
        print('file exists: ', chk_fumia_with_ccle_data_cell_mean)
        return 

    pert_res = pd.read_json(chk_pert_res).transpose()
    query_drugs = pd.read_csv(dataset_query_drugs)
    ccle_dose_resp = pd.read_csv(dataset_ccle_dose_resp)
    ccle_crc = pd.read_csv(dataset_ccle_crc_info)
    cells_df = ccle_crc[['Cell line primary name']].drop_duplicates()
    ccle_dose_resp_sel = ccle_dose_resp.merge(cells_df, how='inner', left_on='Primary Cell Line Name', right_on='Cell line primary name' )
    # ccle_dose_resp_sel = ccle_dose_resp
    ccle_drug_chemblid = query_drugs[['compound_name', 'compound_chembl_id']].drop_duplicates()

    # with compound_chembl_id 
    ccle_dose_resp_sel_2 = ccle_dose_resp_sel.merge(ccle_drug_chemblid, how='inner', left_on='Compound', right_on='compound_name')
    ccle_dose_resp_sel_3 = ccle_dose_resp_sel_2.merge(pert_res, how='inner', left_on='compound_chembl_id', right_index=True)
    ccle_dose_resp_sel_3.to_csv(chk_fumia_with_ccle_data, index=False)
    ccle_dose_resp_sel_3.groupby('Compound').mean().to_csv(chk_fumia_with_ccle_data_cell_mean)


def test_graph_generator(with_small, force):
    
    from sklearn import linear_model
    from sklearn.metrics import r2_score

    df0 = pd.read_csv(chk_fumia_with_ccle_data_cell_mean)            

    # X = df0[['A','P','Q']]    
    # y = df0['ActArea']    
    # reg = linear_model.LinearRegression(fit_intercept=False)    
    # reg.fit (X, y)
    # yhat = reg.predict(X)    
    # r2 = r2_score(y, yhat)
    # print('r2=', r2)

    # df0['yz'] = yz
    # df0['yhat'] = yhat
    # df0.to_csv(output_fumia_with_ccle_data_cell_mean_yhat,index=False)
    
    df0.plot(x='ActArea', y='P', kind='scatter');
    plt.savefig('figure-actarea-vs-p.png'); plt.close()
    
    df0.plot(x='ActArea', y='Q', kind='scatter'); 
    plt.savefig('figure-actarea-vs-q.png'); plt.close()
    
    df0.plot(x='ActArea', y='A', kind='scatter'); 
    plt.savefig('figure-actarea-vs-a.png'); plt.close()




dataset_ccle_crc_info = 'dataset-ccle_crc_info.csv'
dataset_query_drugs = '../table_s2/dataset-query-drugs.csv'
dataset_ccle_dose_resp = '../table_s2/download/CCLE_NP24.2009_Drug_data_2015.02.24.csv'
chk_target_table = 'chk-targetting-table.csv'
chk_drug_model_target_map = 'chk-drug-modeltarget-map.csv'
chk_drug_model_target_dict = 'chk-drug-modeltarget.json'
chk_pert_res = 'chk-pert-result.json'
chk_fumia_with_ccle_data = 'chk-a-fumia-with-ccle-data.csv'
chk_fumia_with_ccle_data_cell_mean = 'chk-a-fumia-with-ccle-data-cells-mean.csv'
output_fumia_with_ccle_data_cell_mean_yhat = 'output-a-fumia-with-ccle-data-cells-mean-yhat.csv'

max_num_drugs = 3;

mut_configs = [
    {
        'State_APC': 'off',
        'State_Ras': 'on',
        'State_Smad': 'off',
        'State_PTEN': 'off',
        'State_p53': 'off',
    }
]



