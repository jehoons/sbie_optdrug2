from os.path import dirname,exists,join
from urllib import request 
from pdb import set_trace
import pandas as pd 
import urllib
import gzip, shutil
import pybel
import pytest 
import json
from os.path import dirname,join,exists
from sbie_optdrug2.results import table_s2 
import sqlite3 
import pickle
import numpy as np 
from tqdm import tqdm
from multiprocessing import Pool


def test_step_1(with_small, force):
    
    if exists(chk_integ_dataset) and force==False:
        print('output already exists:')
        print(chk_integ_dataset)
        return

    conn = sqlite3.connect(table_s2.dataset_chembl_db)
    cur = conn.cursor()

    sql1 = """SELECT 
        MOLECULE_DICTIONARY.CHEMBL_ID 
            as compound_chembl_id,
        COMPOUND_STRUCTURES.CANONICAL_SMILES,
        compound_records.COMPOUND_KEY, 
        compound_records.COMPOUND_NAME,
        target_dictionary.chembl_id
            as target_chembl_id, 
        target_dictionary.pref_name
            as target_name
        FROM compound_records 
        INNER JOIN drug_mechanism 
            on compound_records.RECORD_ID = drug_mechanism.RECORD_ID 
        LEFT JOIN target_dictionary 
            on target_dictionary.TID = drug_mechanism.TID 
        LEFT JOIN MOLECULE_DICTIONARY 
            on compound_records.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO 
        LEFT JOIN COMPOUND_STRUCTURES 
            on COMPOUND_STRUCTURES.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO"""

    if with_small: 
        sql1 += ' LIMIT 100'

    allWithTarget = pd.read_sql_query(sql1, conn)    
    ccleComp = pd.read_csv(table_s2.dataset_drugtarget_info)

    set1 = set(allWithTarget['compound_chembl_id'].unique())
    set2 = set(pd.read_csv(table_s2.dataset_drugtarget_info)['compound_chembl_id'].unique())

    sql2 = """SELECT DISTINCT
        MOLECULE_DICTIONARY.CHEMBL_ID 
            as compound_chembl_id, 
        COMPOUND_STRUCTURES.CANONICAL_SMILES
        FROM MOLECULE_DICTIONARY 
        LEFT JOIN compound_records 
            on compound_records.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO 
        LEFT JOIN COMPOUND_STRUCTURES 
            on COMPOUND_STRUCTURES.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO"""
    
    sql2 += ' WHERE MOLECULE_DICTIONARY.CHEMBL_ID IN (%s)' % \
            ",".join(["'%s'" % c for c in (set2 - set1)])

    df2 = pd.read_sql_query(sql2, conn)    
    df3 = pd.concat([allWithTarget, df2], ignore_index=True)
    df3['source'] = 'CHEMBL'

    df4 = ccleComp[[
        'compound_name',
        'compound_chembl_id',
        'target_name',
        'target_chembl_id']
        ]

    df4['source'] = 'USER'
    df5 = pd.concat([df3, df4], ignore_index=True)
    df5['canonical_smiles'] = df5['canonical_smiles'].fillna('')
    df5['target_chembl_id'] = df5['target_chembl_id'].fillna('')

    idx_drop = (df5['canonical_smiles'] == '') & (
            df5['target_chembl_id'] == '')

    df5[ ~idx_drop ]

    df5.to_csv(chk_integ_dataset, index=False)


def test_step_2(with_small, force):
    
    if exists(chk_similarity_mat) and force==False:
        print('output already exists:')
        print(chk_similarity_mat)
        return 

    dfIntegDB = pd.read_csv(chk_integ_dataset)
    dfIntegDBGrp = dfIntegDB.groupby('compound_chembl_id').groups

    idSmiMap = {}; chembl_ids = [] 
    for g in dfIntegDBGrp: 
        smiles = dfIntegDB.loc[dfIntegDBGrp[g][0], 'canonical_smiles']
        if type(smiles) == float: 
            smiles = None 

        idSmiMap[g] = smiles; chembl_ids.append(g)


    print('Prepare input data for similarity matrix ...')
    inputs = []
    sim_mat = np.zeros([len(chembl_ids), len(chembl_ids)])
    for i,s1 in enumerate(tqdm(chembl_ids)): 
        # progressbar
        for j,s2 in enumerate(chembl_ids):            
            inputs.append([i,j,idSmiMap[s1],idSmiMap[s2]])

    ''' Now, we run parallel processing with multiprocessing ''' 
    
    print('Compute similarity matrix ...')
    pool = Pool(n_cpus); results = []
    for x in tqdm(pool.imap(simmat_worker, inputs), total=len(inputs)):
        results.append(x)

    simmat = np.reshape(results, [len(chembl_ids),len(chembl_ids)])    

    with open(chk_similarity_mat, 'wb') as fobj:
        pickle.dump(simmat, fobj)

    pd.DataFrame(chembl_ids, columns=['compound_chembl_id']
                ).to_csv(chk_similarity_mat_label,index=False)


def simmat_worker(inp):
        
    ''' calculate similarity matrix '''

    i = inp[0]
    j = inp[1]
    s1= inp[2]
    s2= inp[3]

    if i == j: 
        return 1.0
    else:
        if (s1 == None) or (s2 == None) :
            return 0.0
        else: 
            m1 = pybel.readstring("smiles", s1)
            m2 = pybel.readstring("smiles", s2)
            return m1.calcfp() | m2.calcfp()


def test_step_3(with_small, force):                
    
    ''' discover alternative targets ''' 

    if exists(table_s2.output_a_search_res) and force==False:
        print('output already exists:')
        print(table_s2.output_a_search_res)
        return

    search_res = {} 
    for similarity_thr in similarity_list:        
        keystr = 'similarity>=%.02f' % similarity_thr
        search_res[keystr] = resembl_worker(similarity_thr)

    with open(table_s2.output_a_search_res, 'w') as fobj: 
        json.dump(search_res, fobj, indent=4)


def resembl_worker(similarity_thr):

    labels = pd.read_csv(
        chk_similarity_mat_label)['compound_chembl_id'].tolist()
    user_drugset = pd.read_csv(dataset_drugtarget_info) 
    query_list = user_drugset['compound_chembl_id'].unique().tolist()    
    int_db = pd.read_csv(chk_integ_dataset)

    with open(chk_similarity_mat, 'rb') as fobj:
        simmat = pickle.load(fobj)

    f0 = {'compound_name': lambda x: x.unique()[0] \
        if len(x.unique())==1 else 'multiple value error'}

    chemblid_name_map = user_drugset.groupby(
        'compound_chembl_id').agg(f0).to_dict()['compound_name']

    search_res = {} 
    for q0 in tqdm(query_list): # for each query compound 
        # 1. check similar compound
        ix = labels.index(q0)
        selected = simmat[ix, ] >= similarity_thr 
        neighbor = []        
        neighbor_score = {} 
        for i,sel in enumerate(selected):
            if sel: 
                neighbor.append( labels[i] )

        # 2. check targets of neighbor(similar compound)
        neighbor_info = {}; 
        connectivity = {}; 
        connectivity_w = {}
        for nei in neighbor:
            target_neighbor_similarity = simmat[ix, labels.index(nei)]
            subset_idx = (int_db['compound_chembl_id'] == nei) & \
                (~int_db['target_chembl_id'].isnull())
            target_info = int_db[subset_idx][
                ['target_name','target_chembl_id','source']].to_dict(orient='records')
            neighbor_info[nei] = {
                'similarity': target_neighbor_similarity, 
                'target_info': target_info
                }
            for tar in target_info:
                chemblid = tar['target_chembl_id']
                if chemblid not in connectivity:
                    connectivity[chemblid] = 1
                else: 
                    connectivity[chemblid] +=1

                if chemblid not in connectivity_w:
                    connectivity_w[chemblid] = target_neighbor_similarity
                else: 
                    connectivity_w[chemblid] +=target_neighbor_similarity

        search_res[q0] = {
                'neighbor_info': neighbor_info, 
                'connectivity': connectivity, 
                'connectivity_w': connectivity_w,
                'compound_name': chemblid_name_map[q0]
                }
                
    return search_res


def test_step_4(with_small, force):
    
    ''' postprocess ''' 

    with open(table_s2.output_a_search_res, encoding='utf-8') as f: 
        search_res = json.loads(f.read())

    mykey = 'similarity>=0.35'

    res = search_res[mykey]

    dforig = pd.read_csv(table_s2.dataset_drugtarget_info)
    ccle = dforig['target_chembl_id'].unique().tolist()

    df0 = pd.read_csv(table_s2.dataset_model_node_info)    
    fumia_nodes = df0['chembl_id'].dropna().unique().tolist()
    
    targets = []
    for comp in res: 
        for tar in res[comp]['connectivity_w']:
            tar_score = res[comp]['connectivity_w'][tar]
            if tar_score >= 0.01:
                targets.append(tar)

    set_fumia_nodes = set(fumia_nodes)

    dforig['targetable'] = 0 
    for i in dforig.index: 
        if dforig.loc[i, 'target_chembl_id'] in set_fumia_nodes:
            dforig.loc[i, 'targetable'] = 1

    df_ccle_extended = pd.DataFrame([], 
        columns=[
            'compound_chembl_id',
            'target_chembl_id',
            'connectivity_w',
            'targetable'
            ]
        )

    i = 0
    for compound_chembl_id in res: 
        for target_chembl_id in res[compound_chembl_id]['connectivity_w']: 
            df_ccle_extended.loc[i, 'compound_chembl_id'] = compound_chembl_id
            df_ccle_extended.loc[i, 'target_chembl_id'] = target_chembl_id
            df_ccle_extended.loc[i, 'connectivity_w'] = \
                res[compound_chembl_id]['connectivity_w'][target_chembl_id]
            if df_ccle_extended.loc[i, 'target_chembl_id'] in set_fumia_nodes:
                df_ccle_extended.loc[i, 'targetable'] = 1
            else: 
                df_ccle_extended.loc[i, 'targetable'] = 0

            i += 1

    fcn_set = {
        'targetable': {
            'count': 'count',
            'coverage%': lambda g: g.sum()/g.count()*100.0
            }
        }
    
    dforig.groupby('compound_chembl_id').agg(
            fcn_set).to_csv(chk_drug_coverage)
    
    df_ccle_extended.groupby('compound_chembl_id').agg(
            fcn_set).to_csv(chk_drug_coverage_ext)


chk_integ_dataset = join(table_s2.basedir, 'chk-integ-dataset.csv')

chk_similarity_mat = join(table_s2.basedir, 'chk-similarity-mat.pkl')

chk_similarity_mat_label = join(table_s2.basedir, 'chk-similarity-mat-label.csv')

# coverage for original information
chk_drug_coverage = join(table_s2.basedir, 'chk-drug-coverage.csv')

# coverage for extended information
chk_drug_coverage_ext = join(table_s2.basedir, 'chk-drug-coverage-ext.csv')

# min_similarity = 0.45
similarity_list = [0.30,0.35,0.5,1.0]

n_cpus = 100

