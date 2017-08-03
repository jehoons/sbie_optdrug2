from os.path import dirname,exists,join
import urllib
from urllib import request 
import gzip, shutil

basedir = dirname(__file__)

dataset_chembl_db = join(basedir, 
    'download/chembl_23/chembl_23_sqlite/chembl_23.db')

dataset_chembl_sdf = join(basedir,
    'download/chembl_23.sdf')

dataset_chembl_uniprot = join(basedir,
    'download/chembl_uniprot_mapping.txt')

dataset_ccle_treatment = join(basedir,
    'download/CCLE_NP24.2009_Drug_data_2015.02.24.csv')

# input query drugs
dataset_drugtarget_info = join(basedir, 
    'dataset-query-drugs.csv')

# input model node information
dataset_model_node_info = join(basedir, 
    'dataset-fumia-node-info-update-2.csv')

output_a_search_res = join(basedir, 
    'output-a-alternative-targets.json')

