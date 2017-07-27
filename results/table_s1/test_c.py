import pybel, pandas as pd, pytest, json, sqlite3, pickle, numpy as np
from tqdm import tqdm 
from ipdb import set_trace
from os.path import dirname,join,exists
from sbie_optdrug2.results import table_s1


def test_c1(force, with_small):

    ''' summary CCLE compounds '''

    if exists(table_s1.output_c1) and force==False:
        return 

    df = pd.read_csv(table_s1.dataset_ccle_treatment)
    df[['Compound','Target']].drop_duplicates().to_csv(
        table_s1.output_c1, index=False)


def test_c2(force, with_small):

    ''' Discover alternative targets based on compound similarity '''

    df = pd.read_csv(table_s1.updated_c1)    
    query_list = df['compound_chembl_id'].unique().tolist()
    res = find_resemble_compounds(
        query_list=query_list, min_similarity=0.2)    
    with open(table_s1.output_c2, 'w') as f: 
        json.dump(res, f, indent=4)


def find_resemble_compounds(query_list=[], min_similarity=0.5): 

    with open(table_s1.output_b2_drug_simil, 'rb') as f: 
        data = pickle.load(f)
    
    dataint = pd.read_csv(table_s1.output_b1_data_integration)    
    labelinfo = pd.read_csv(table_s1.output_b3_drug_simil_labels)     
    labels = labelinfo['chembl_id_x'].tolist()

    output = {}
    for q0 in query_list: 
        if q0 not in labels: 
            print('query(%s) not found in the matrix' % q0)
            continue 

        ix = labels.index(q0)        
        selected = data[ix, ] > min_similarity 
        neighbor = []        
        
        for i,sel in enumerate(selected):
            if sel: 
                neighbor.append( labels[i] )

        chembl_grp = dataint.groupby('chembl_id_x').groups

        comp_counter = {}        
        neighbor_dict = {}
        
        for nei in neighbor:             
            jx = labels.index(nei)
            neighbor_dict[nei] = data[ix, jx]
            ix2 = chembl_grp[nei]
            for target in dataint.loc[ix2, 't_chembl_id'].tolist():                
                if pd.isnull(target): 
                    continue 
                if target not in comp_counter:
                    comp_counter[target] = 1
                else: 
                    comp_counter[target] += 1

        output[q0] = {
            'inferred_targets': comp_counter, 
            'similar_compounds': neighbor_dict
            }

    return output

    # to dataframe? 
    