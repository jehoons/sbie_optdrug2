from os.path import dirname,exists,join
import urllib
from urllib import request 
import gzip, shutil


basedir = dirname(__file__)

# dataset
dataset_chembl_db = join(basedir, 
    'download/chembl_23/chembl_23_sqlite/chembl_23.db')

dataset_chembl_sdf = join(basedir, 
    'download/chembl_23.sdf')

dataset_chembl_uniprot = join(basedir,
    'download/chembl_uniprot_mapping.txt')

dataset_ccle_treatment = join(basedir,
    'download/CCLE_NP24.2009_Drug_data_2015.02.24.csv')

# output files
output_a1 = join(basedir, 
    'output-a1-filt-1st-with-target-info.csv')

output_a2 = join(basedir, 
    'output-a2-uniprot-names.csv')

# uniprot.org를 통해서 유전자 이름으로 변환할 수가 있다. 
handwork_a2 = join(basedir, 
    'handwork-a2-gene-names.csv')

output_a3 = join(basedir, 
    'output-a3-compounds.csv')

output_b1_data_integration = join(basedir, 
    'output-b1-data-integration.csv')

output_b2_drug_simil = join(basedir, 
    'output-b2-drug-similarity.pkl')

output_b2_drug_simil_csv = join(basedir, 
    'output-b2-drug-similarity.csv')

output_b3_drug_simil_labels = join(basedir, 
    'output-b3-drug-simil-labels.csv')

output_b4 = join(basedir, 
    'output-b4-comp-target-dict.csv')

output_c1 = join(basedir, 
    'output-c1-ccle-drug-target-pair.csv')

handwork_c1 = join(basedir, 
    'handwork-c1-ccle-drug-target-pair.csv')

output_c2 = join(basedir, 
    'output-c2-inferred-targets.json')

output_c2_stats = join(basedir, 
    'output-c2-inferred-targets-stats.csv')

output_c3 = join(basedir, 
    'output-c3')
