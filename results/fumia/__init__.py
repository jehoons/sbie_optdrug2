import json,itertools,pytest
from os.path import exists,dirname,join
from numpy.random import random
from boolean3_addon import attr_cy, to_logic
# from pdb import set_trace
# from boolean3_addon import to_logic 
# import pandas as pd
# from tqdm import tqdm

# from sbie_optdrug.result_season2 import tab_s1
# from sbie_optdrug.result_season2.tab_s1 import fumia_phenotype
# from fumia_phenotype import *

# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt

__all__ = [] 

def modeltext():

    datast_fumia = join(dirname(__file__), 'dataset-fumia-model-ws.txt')    
    with open(datast_fumia, 'r') as fobj:
        lines = fobj.readlines()
        lines2 = [] 
        for lin in lines: 
            lin = lin.strip()
            if lin[0] == '#': 
                continue 
            right = lin.split('=')[1].strip()
            if right == 'input':            
                lines2.append( lin.split('=')[0].strip() + '=' + 'False') 
            else: 
                lines2.append(lin)

    return "\n".join(lines2)


def get_input_combinations():    

    inputs = [ {
    'State_Mutagen' : State_Mutagen, 
    'State_GFs': State_GFs,
    'State_Nutrients': State_Nutrients,
    'State_TNFalpha': State_TNFalpha,
    'State_Hypoxia': State_Hypoxia,
    'State_Gli': False,
    } for State_Mutagen,State_GFs,State_Nutrients,State_TNFalpha, \
        State_Hypoxia in itertools.product([False,True], repeat=5)]

    return inputs


def to_bits(input_data):    
    res = ['1' if input_data[state] else '0' for state in ['State_Mutagen','State_GFs','State_Nutrients','State_TNFalpha','State_Hypoxia','State_Gli']]
    
    return "".join(res)


def test_to_bits():
    inp = get_input_combinations()
    res = to_bits(inp[0])
    
    print(res)


def get_mutations_config():
   
    # free
    mutations_config_list = [ 
        {
            'State_APC': None,
            'State_Ras': None,
            'State_Smad': None,
            'State_PTEN': None,
            'State_p53': None,
        },

        {
            'State_APC': 'off',
            'State_Ras': None,
            'State_Smad': None, 
            'State_PTEN': None, 
            'State_p53': None, 
        }, 

        {
            'State_APC': 'off',
            'State_Ras': 'on',
            'State_Smad': None, 
            'State_PTEN': None, 
            'State_p53': None, 
        }, 

        {
            'State_APC': 'off',
            'State_Ras': 'on',
            'State_Smad': 'off',
            'State_PTEN': None, 
            'State_p53': None, 
        }, 
        
        {
            'State_APC': 'off',
            'State_Ras': 'on',
            'State_Smad': 'off',
            'State_PTEN': 'off', 
            'State_p53': None, 
        },

        {
            'State_APC': 'off',
            'State_Ras': 'on',
            'State_Smad': 'off',
            'State_PTEN': 'off', 
            'State_p53': 'off',
        }
    ]
    return mutations_config_list


def hello():

    if not exists('fumia_engine.pyx'):
        attr_cy.build(modeltext(),
            pyx='fumia_engine.pyx',
            weighted_sum=False
            )

    import pyximport; pyximport.install()
    import fumia_engine

    res = attr_cy.parallel(
        fumia_engine, steps=80, samples=100, ncpus=100, 
        repeats=1000, on_states=[], off_states=[]
        ) 

    with open(output_hello, 'w') as f: 
        json.dump(res, f, indent=4)


def test_hello():
    
    hello()


basedir = dirname(__file__)
   
output_hello = 'chk-hello-fumia-attr.json'

dataset_fumia_node_info = join(basedir, 'dataset-fumia-node-info-update-2.csv')