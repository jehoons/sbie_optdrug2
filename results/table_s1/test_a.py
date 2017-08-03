from pdb import set_trace
from boolean3_addon import to_logic 
from os.path import exists,dirname
from sbie_optdrug.result_season2 import tab_s1
import pandas as pd 

def test_1():

    with open(file_a1, 'r') as fobj: 
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

    with open(file_a2, 'w') as fobj:
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

    with open(file_a2, 'r') as fobj:
        lines = fobj.readlines()

    res = to_logic.build("".join(lines), short=False)

    with open(file_a3, 'w') as fobj2:
        fobj2.write(res)

    # set_trace()

def test_2():
    with open(file_a3, 'r') as fobj: 
        lines = fobj.readlines()

    state_list = [] 
    for l0 in lines: 
        l0 = l0.strip()
        if l0 == '' : continue 
        if l0[0] == '#' : continue

        left = l0.split('=')[0].replace('*','')
        left = left.replace(' ','')
        # left = left.split('_')[1]
        # print(left)
        state_list.append(left)

    state_list = [x for x in set(state_list)]

    from tqdm import tqdm
    import json
    data = json.load(open('a/hgnc_complete_set.json','r'))
    docs = data['response']['docs']

    res = {} 

    for state in tqdm(state_list):
        state2 = state.replace('State_', '')
        # print(state2)
        guess_list = []
        state_parts = state2.split('_')
        for doc in (docs): 
            docstr = repr(doc)
            # print (docstr)
            if docstr.find(state2) >= 0: 
                # print (docstr)
                guess = doc['symbol']
                guess_list.append(guess)

            # for part in state_parts: 
            #     if docstr.find(part) >= 0:                 
            #         guess = doc['symbol']
            #         guess_list.append(guess)

        res[state] = guess_list

        # break; 


    json.dump(res, open(file_a5, 'w'), indent=4)
    # df0.to_csv(file_a5)
    # set_trace()


def test_3():
    # import json
    # data = json.load(open('a/hgnc_complete_set.json','r'))

    # (Pdb++) data['response']['docs'][0]
    # {'omim_id': ['138670'], 'locus_group': 'protein-coding gene', 
    # 'rgd_id': ['RGD:69417'], 'location': '19q13.43', 'ccds_id': ['CCDS12976'], 
    # 'vega_id': 'OTTHUMG00000183507', 'gene_family': ['Immunoglobulin like domain containing'],
    # 'ensembl_gene_id': 'ENSG00000121410', 'pubmed_id': [2591067], 'symbol': 'A1BG', 
    # 'status': 'Approved', 'mgd_id': ['MGI:2152878'], 'gene_family_id': [594], 
    # 'date_approved_reserved': '1989-06-30', 'locus_type': 'gene with protein product', 
    # 'ucsc_id': 'uc002qsd.5', 'location_sortable': '19q13.43', 'name': 'alpha-1-B glycoprotein',
    # 'entrez_id': '1', 'merops': 'I43.950', 'refseq_accession': ['NM_130786'], 
    # 'uuid': '82e66fd2-089d-4f2e-8610-320c9f7baf5b', '_version_': 1564445527844585472, 'uniprot_ids': ['P04217'], 
    # 'date_modified': '2015-07-13', 'cosmic': 'A1BG', 'hgnc_id': 'HGNC:5'}
    # (Pdb++) data['response']['docs'][0]

    # docs = data['response']['docs']


    # set_trace()

    pass 



# input 
file_a1 = dirname(tab_s1.__file__) + '/a/a1-fumia-model-curated-from-suppl.txt'

# output 
file_a2 = dirname(tab_s1.__file__) + '/a/a2-fumia-model-processed-weighted-sum.txt'
file_a3 = dirname(tab_s1.__file__) + '/a/a3-fumia-model-processed-logic.txt'

file_a5 = dirname(tab_s1.__file__) + '/a/a5-fumia-model-nodes.json'