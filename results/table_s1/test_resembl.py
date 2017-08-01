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
from sbie_optdrug2.results import table_s1 
import sqlite3 
import pickle
import numpy as np 
from tqdm import tqdm
from multiprocessing import Pool

# sql = """SELECT 
#     MOLECULE_DICTIONARY.CHEMBL_ID,
#     compound_records.RECORD_ID,  
#     compound_records.MOLREGNO,  
#     compound_records.COMPOUND_KEY,  
#     compound_records.COMPOUND_NAME, 
#     drug_mechanism.TID, 
#     COMPOUND_STRUCTURES.CANONICAL_SMILES,
#     -- COMPOUND_STRUCTURES.MOLFILE,
#     target_dictionary.pref_name as t_name, 
#     target_dictionary.chembl_id as t_chembl_id 
#     FROM compound_records 
#     INNER JOIN drug_mechanism on compound_records.RECORD_ID = drug_mechanism.RECORD_ID 
#     LEFT JOIN target_dictionary on target_dictionary.TID = drug_mechanism.TID 
#     LEFT JOIN MOLECULE_DICTIONARY on compound_records.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO 
#     LEFT JOIN COMPOUND_STRUCTURES on COMPOUND_STRUCTURES.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO 
# """

def test_step_1(with_small, force):
    
    if exists(chk_integ_dataset) and force==False:
        return
        
    conn = sqlite3.connect(dataset_chembl_db)
    cur = conn.cursor()

    sql1 = """SELECT 
        MOLECULE_DICTIONARY.CHEMBL_ID 
            as compound_chembl_id,
        COMPOUND_STRUCTURES.CANONICAL_SMILES,
        target_dictionary.chembl_id
            as target_chembl_id 
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
        sql1 += ' LIMIT 10'

    allWithTarget = pd.read_sql_query(sql1, conn)    
    ccleComp = pd.read_csv(dataset_drugtarget_info)

    set1 = set(allWithTarget['compound_chembl_id'].unique())
    set2 = set(pd.read_csv(dataset_drugtarget_info)['compound_chembl_id'].unique())

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

    df4 = ccleComp[['compound_chembl_id', 'target_chembl_id']]
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
        return 

    dfIntegDB = pd.read_csv(chk_integ_dataset)
    dfIntegDBGrp = dfIntegDB.groupby('compound_chembl_id').groups

    idSmiMap = {}; chembl_ids = [] 
    for g in dfIntegDBGrp: 
        smiles = dfIntegDB.loc[dfIntegDBGrp[g][0], 'canonical_smiles']
        if type(smiles) == float: 
            smiles = None 

        idSmiMap[g] = smiles; chembl_ids.append(g)

    inputs = []

    sim_mat = np.zeros([len(chembl_ids), len(chembl_ids)])
    for i,s1 in enumerate(tqdm(chembl_ids)):
        for j,s2 in enumerate(chembl_ids):            
            inputs.append([i,j,idSmiMap[s1],idSmiMap[s2]])

    ''' Now, we run parallel processing with multiprocessing ''' 
    
    pool = Pool(n_cpus); results = []

    for x in tqdm(pool.imap(simmat_worker, inputs), total=len(inputs)):
        results.append(x)

    simmat = np.reshape(results, [len(chembl_ids),len(chembl_ids)])    

    with open(chk_similarity_mat, 'wb') as fobj:
        pickle.dump(simmat, fobj)

    pd.DataFrame(chembl_ids, columns=['compound_chembl_id']
                ).to_csv(chk_similarity_mat_label,index=False)


def simmat_worker(inp):
        
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

    if exists(output_search_res) and force==False:
        return

    labels = pd.read_csv(chk_similarity_mat_label)['compound_chembl_id'].tolist()
    user_drugset = pd.read_csv(dataset_drugtarget_info) 
    query_list = user_drugset['compound_chembl_id'].unique().tolist()    
    int_db = pd.read_csv(chk_integ_dataset)

    with open(chk_similarity_mat, 'rb') as fobj:
        simmat = pickle.load(fobj)

    search_res = {} 
    for q0 in tqdm(query_list): 
        ix = labels.index(q0)
        selected = simmat[ix, ] > min_similarity 
        neighbor = []        
        neighbor_score = {} 

        for i,sel in enumerate(selected):
            if sel: 
                neighbor.append( labels[i] )

        neighbor_info = {}
        for nei in neighbor:
            subset_idx = (int_db['compound_chembl_id'] == nei) & (
                    ~int_db['target_chembl_id'].isnull())

            subset = int_db[subset_idx]

            target_info = subset[['target_chembl_id','source']].to_dict(
                        orient='records')

            neighbor_info[nei] = {
                'similarity': simmat[ix, labels.index(nei)], 
                'target': target_info
                }

        connectivity = {} 
        for tar in target_info:
            chemblid = tar['target_chembl_id']
            if chemblid not in connectivity:
                connectivity[chemblid] = 1 
            else: 
                connectivity[chemblid] +=1 

        search_res[q0] = {
                'neighbor_info': neighbor_info, 
                'target_info': target_info, 
                'connectivity': connectivity
                }
        
    with open(output_search_res, 'w') as fobj: 
        json.dump(search_res, fobj, indent=1)


basedir = dirname(__file__)

# background dataset 
dataset_chembl_db = table_s1.dataset_chembl_db

# dataset for input 
dataset_drugtarget_info = join(basedir, 
        'dataset-query-drugs.csv')

# dataset for output 
output_search_res = join(basedir, 
        'output-a-alternative-targets.json')

# datafiles for check 
chk_integ_dataset = join(basedir, 
        'chk-integ-dataset.csv')

chk_similarity_mat = join(basedir, 
        'chk-similarity-mat.pkl')

chk_similarity_mat_label = join(basedir, 
        'chk-similarity-mat-label.csv')

min_similarity = 0.3

n_cpus = 100

