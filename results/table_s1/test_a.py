from pdb import set_trace
from boolean3_addon import to_logic 
from os.path import exists,dirname
from sbie_optdrug2.results import table_s1
import pandas as pd 

def test_1():

    with open(dataset_fumianet, 'r') as fobj: 
        lines = fobj.readlines()

    lines2 = [] 
    for l in lines: 
        l = l.strip() 
        if l == '': continue 
        if l[0] == '#': continue     
        l = l.replace(' ', '')
        lines2.append(l)

    states = [] 
    for l in lines2: 
        words = l.split('(t+1)=')
        states.append(words[0])

    repdict = {}
    for st in states: 
        st2 = st 
        st2 = st2.replace('α','alpha')
        st2 = st2.replace('β','beta')
        st2 = st2.replace('κ','kappa')
        st2 = st2.replace('-','_')
        st2 = st2.replace('−','_')
        st2 = st2.replace('/','_')
        st2 = st2.replace('/','_')
        st2 = st2.replace('/','_')
        repdict[st.strip()] = st2.strip()

    repdict['sgn['] = 'sign('
    repdict[']'] = ')'
    repdict['(t+1)'] = ''
    repdict['(t)'] = ''
    repdict[';'] = ''

    lines3 = [] 
    for thisline in lines2: 
        # print('before:', thisline)
        for rep in repdict:
            thisline = thisline.replace(rep, repdict[rep])
        thisline = thisline.replace('σ','State_')
        thisline = thisline.replace('−', '-')
        thisline = thisline.replace('/', '_')
        lines3.append(thisline)

    with open(output_a1, 'w') as fobj:
        for lin in lines3:
            words = lin.split('=')
            if words[1].strip() == 'input':
                fobj.write(lin+'\n')
            else:
                fobj.write(words[0]+' = '+'Random'+'\n')

        for lin in lines3:
            words = lin.split('=')
            if words[1].strip() == 'input':
                fobj.write(words[0]+' = '+words[0]+'\n')
            else:
                fobj.write(words[0]+' *= '+words[1]+'\n')

    with open(output_a1, 'r') as fobj:
        lines = fobj.readlines()

    res = to_logic.build("".join(lines), short=False)

    with open(output_a2, 'w') as fobj2:
        fobj2.write(res)

# input 
dataset_fumianet = dirname(table_s1.__file__) + '/dataset-fumia-curated-rev1.txt'

# output 
output_a1 = dirname(table_s1.__file__) + '/output-a1-fumia-model-processed-weighted-sum.txt'
output_a2 = dirname(table_s1.__file__) + '/output-a2-fumia-model-processed-logic.txt'

