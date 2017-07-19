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

    ''' summary CCLE compounds '''

    if exists(table_s1.output_c1) and force==False:
        return 

    df = pd.read_csv(table_s1.dataset_ccle_treatment)
    df[['Compound','Target']].drop_duplicates().to_csv(
        table_s1.output_c1, index=False)
