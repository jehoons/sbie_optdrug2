from pdb import set_trace
import pybel
import pandas as pd 
import pytest 
import json
from os.path import dirname,join,exists
from sbie_optdrug2.results import table_s1 
import sqlite3 
import pickle 

def test_a(with_small, force):

    conn = sqlite3.connect(table_s1.dataset_chembl_db)
    cur = conn.cursor()

    sql = """SELECT 
        MOLECULE_DICTIONARY.CHEMBL_ID,
        compound_records.RECORD_ID,  
        compound_records.MOLREGNO,  
        compound_records.COMPOUND_KEY,  
        compound_records.COMPOUND_NAME, 
        drug_mechanism.TID, 
        COMPOUND_STRUCTURES.CANONICAL_SMILES,
        -- COMPOUND_STRUCTURES.MOLFILE,
        target_dictionary.pref_name as t_name, 
        target_dictionary.chembl_id as t_chembl_id 
        FROM compound_records 
        INNER JOIN drug_mechanism on compound_records.RECORD_ID = drug_mechanism.RECORD_ID 
        LEFT JOIN target_dictionary on target_dictionary.TID = drug_mechanism.TID 
        LEFT JOIN MOLECULE_DICTIONARY on compound_records.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO 
        LEFT JOIN COMPOUND_STRUCTURES on COMPOUND_STRUCTURES.MOLREGNO = MOLECULE_DICTIONARY.MOLREGNO 
        -- LIMIT 10 
    """
    
    df = pd.read_sql_query(sql, conn)

    df_uniprot = pd.read_csv(table_s1.dataset_chembl_uniprot, sep='\t', skiprows=1,
        names=['uniprot_id', 'chembl_id', 'target_desc', 'protein_type'])
    
    df_integrated = pd.merge(df, df_uniprot, how='left', left_on='t_chembl_id', right_on='chembl_id')
    
    df_integrated.to_csv(table_s1.output_a1)
    
    df_integrated.groupby('uniprot_id').count()[[]].to_csv(table_s1.output_a2)
    
    df_integrated.groupby('chembl_id_x').count()[[]].to_csv(table_s1.output_a3)
