from pdb import set_trace
import pybel
import pandas as pd 
import pytest 
import json
from os.path import dirname,join,exists
from sbie_optdrug2.results import table_s1 
import sqlite3 
import pickle 
import numpy as np
from tqdm import tqdm 

def test_b1(force, with_small):
    
    ''' extract moles from sdf file '''
    df0 = pd.read_csv(table_s1.output_a1)    

    df_uniprot_map = pd.read_csv(table_s1.dataset_a3_names)

    df0 = pd.merge(df0, df_uniprot_map, how='left', left_on='uniprot_id', right_on='From')

    df0.to_csv(table_s1.output_b1_data_integration) 

    df_unique = df0[['chembl_id_x', 'canonical_smiles']].drop_duplicates().dropna(how='any')    
    
    smiles = df_unique['canonical_smiles'].tolist()

    if with_small: 
        smiles = smiles[0:20]
    
    sim_mat = np.zeros([len(smiles), len(smiles)])
    
    for i,s1 in enumerate(tqdm(smiles)):
        for j,s2 in enumerate(smiles):
            if sim_mat[j,i] == 0.0: 
                m1 = pybel.readstring("smiles", smiles[i])
                m2 = pybel.readstring("smiles", smiles[j])
                sim_mat[i,j] = m1.calcfp() | m2.calcfp()
            else: 
                sim_mat[i,j] = sim_mat[j,i]
     
    with open(table_s1.output_b2_drug_simil, 'wb') as fobj:
        pickle.dump(sim_mat, fobj)

    pd.DataFrame(sim_mat).to_csv('output-similarity-matrix.csv')

    df_unique.to_csv(table_s1.output_b3_drug_simil_labels)

    comp_target_dict = df0.groupby('chembl_id_x')['To'].apply(lambda x: ",".join([str(q) for q in x]))
    comp_target_dict = comp_target_dict.to_frame()

    comp_target_dict2 = df0.groupby('chembl_id_x')['compound_name'].apply(lambda x: x.unique()[0])
    comp_target_dict['compound_name'] = comp_target_dict2

    comp_target_dict[['compound_name','To']].to_csv(table_s1.output_b4)



