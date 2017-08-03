import os
import sys
from pdb import set_trace
import pandas as pd 
from dataset_ccle import load_treatment, load_cellinfo

df_treat = load_treatment()
df_cellinfo = load_cellinfo()
tissues = set(df_cellinfo['Site Primary'])

# all tissue: 
# {
# 'salivary_gland', 'pleura', 'ovary', 'endometrium', 'large_intestine',
# 'central_nervous_system', 'autonomic_ganglia', 'biliary_tract', 
# 'soft_tissue', 'skin', 'oesophagus', 'stomach', 'prostate', 'breast', 
# 'bone', 'liver', 'pancreas', 'kidney', 'thyroid', 'urinary_tract', 
# 'upper_aerodigestive_tract', 'lung', 'haematopoietic_and_lymphoid_tissue',
# 'small_intestine'
# }

selected_tissue = [
	'large_intestine',
	# 'biliary_tract',
	# 'urinary_tract',
	# 'upper_aerodigestive_tract',
	'small_intestine'
]

for i, tissue in enumerate(selected_tissue):
	if i == 0:
		idx1 = df_cellinfo['Site Primary'] == tissue
	else: 
		idx1 = idx1 | (df_cellinfo['Site Primary'] == tissue)

cells = [x for x in set(df_cellinfo.loc[idx1, 'CCLE name'])]

for i, cell in enumerate(cells):
	if i == 0:
		idx2 = df_treat['CCLE Cell Line Name'] == cell
	else: 
		idx2 = idx2 | (df_treat['CCLE Cell Line Name'] == cell)

# endfor

df_treat_selected = df_treat.loc[idx2]
compounds_selected = [x for x in set(df_treat_selected['Compound'])]
targets_selected = [x for x in set(df_treat_selected['Target'])]

df_treat_selected[['CCLE Cell Line Name', 'Compound', 'Target']].to_csv('output.csv')

# fumia node data 
fumia_node_data_file = 'a/a5-fumia-model-nodes-working.json'
import json 

with open(fumia_node_data_file, 'r') as fin: 
	res = json.load(fin)

data = [] 
for k in res:
	data += res[k]

fumia_genes_extended = [x for x in set(data)]

for tar in set(df_treat_selected['Target']): 
	find_res = -1
	for fum in fumia_genes_extended: 
		find_res = fum.find(tar)
		if find_res == 0: 
			break 
	print (tar, ',', find_res==0)


df_treat_selected['Targetable'] = False 

for x in ['RTK', 'HSP90', 'HDAC', 'ABL', 'FGFR', 'CDK4', 'RAF', 'MDM2']: 
	ix = df_treat_selected['Target'] == x 
	df_treat_selected.loc[ix, 'Targetable'] = True

# 224/535

set_trace()


