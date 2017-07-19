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


def test_c1(force, with_small):
    # ''' extrace ccle drugs '''
    if exists(table_s1.output_c1) and force==False:
        return 

    df = pd.read_csv(table_s1.dataset_ccle_treatment)
    df[['Compound','Target']].drop_duplicates().to_csv(
        table_s1.output_c1, index=False)

# def test_c2(with_small, force):

#     conn = sqlite3.connect(table_s1.dataset_chembl_db)
#     cur = conn.cursor()
#     sql = "SELECT * from MOLECULE_SYNONYMS WHERE SYNONYMS LIKE \'%s\'" 
#     df_ccle_comp = pd.read_csv(table_s1.output_c1)
#     comps = df_ccle_comp['Compound'].tolist()
#     output_dict = {}
#     df1 = pd.DataFrame([])
#     frames = []

#     for c in comps:
#         sql0 = sql % ('%'+c+'%')
#         df = pd.read_sql_query(sql0, conn)

#         syn_list = df['synonyms'].unique().tolist()
#         output_dict[c] = ",".join(syn_list) if len(syn_list)>0 else None

#         if syn_list == []: 
#             df.loc[0] = None 
        
#         df['input_drug'] = c

#         # df['molregno'] = df['molregno'].astype(int)

#         frames.append(df)

#     # set_trace()
#     res = pd.concat(frames, ignore_index=True)
#     res.to_csv('output.csv', index=False)
#     set_trace()