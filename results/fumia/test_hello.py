import json,itertools,pytest
from os.path import exists
from numpy.random import random
from boolean3_addon import attr_cy, to_logic
from pdb import set_trace
from boolean3_addon import to_logic 
from os.path import exists,dirname
import pandas as pd
from tqdm import tqdm

from sbie_optdrug.result_season2 import tab_s1
from sbie_optdrug.result_season2.tab_s1 import fumia_phenotype
from fumia_phenotype import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def preproc_model_file(datast_fumia):

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

    modeltext = "\n".join(lines2)
    return modeltext

    
def test_hello(with_small, force):

    if not exists('fumia_engine.pyx'):
        attr_cy.build(preproc_model_file(datast_fumia),
            pyx='fumia_engine.pyx',
            weighted_sum=False
            )

    import pyximport; pyximport.install()
    import fumia_engine

    res = attr_cy.parallel(
        fumia_engine, steps=80, samples=100, ncpus=100, 
        repeats=1000, on_states=[], off_states=[]
        ) 

    with open(output_attr, 'w') as f: 
        json.dump(res, f, indent=4)

datast_fumia = 'dataset-fumia-model.txt'
output_attr = 'output-a-attractors.json'

