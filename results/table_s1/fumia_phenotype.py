from os.path import dirname,join

def remove_dup(test_list):
    prev = '' 
    output_list = [] 
    for el in test_list:
        if el == prev: 
            continue 
        else: 
            output_list.append(el)
            prev = el 

    return output_list


def check_cyclin(a_state_in_cyc, lbls): 
    global cycs
    rr = [lbls.index(s) for s in cycs]
    res0 = [ a_state_in_cyc[si] for k,si in enumerate(rr)]
    return "".join(res0)


def check_cyclin_seq(att_value, state_key, lbls):
    global cycs
    cyc_att = [state_key[x] for x in att_value]
    cyc_seq = [] 
    for S in cyc_att: 
        res0 = check_cyclin(S,lbls)
        cyc_seq.append("".join(res0))

    rem = remove_dup(cyc_seq)

    return rem


def check_cyclic_phenotype(seq):
    if len(seq) == 5: 
        for i,j in zip(seq, ['0010', '0011', '1011', '1010', '1110']):
            if i != j:
                return 'Q'
        else :
            return 'P'

    elif len(seq) == 8:
        for i,j in zip(seq, ['0011', '1011', '1010', '1110', '0010', '1010', '1110', '0010']):
            if i != j:
                return 'Q'
        else :
            return 'P'

    else: 
        return 'Q'


def attach_phenotype(result):
# from the paper: 
# Considering the effects of mutations, reported on the next subsection, 
# they include the following basic
# cell phenotypes: 
# apoptotic, characterized by active caspases;
# glycolytic, with H1f1 activated under normoxia; 
# immortalized in which hTert is active; 
# migratory, associated to inactivate Ecadherin;
# mutator, corresponding to inactive Atm/Atr proteins in the presence of DNA damage; 
# proliferative, in which cyclins are activated along the cell cycle in the correct sequence; 
# and quiescent, with cyclins inactive or activated in a wrong sequence.
# 위키를 살펴보면, 
# https://en.wikipedia.org/wiki/Cyclin
# D는 항상 켜져 있어야 하고, 
# E, A, B 순서로 켜져야 한다는 것을 알 수 있다.

    simul_result = result['simul_result']
    # input_cond = result['input_condition']

    atts = simul_result['attractors']
    lbls = simul_result['labels']
    state_key = simul_result['state_key']

    for att in atts:            
        att_value = atts[att]['value']
        att_type = atts[att]['type']
        if att_type == 'cyclic':
            seq = check_cyclin_seq(att_value, state_key, lbls)
            result['simul_result']['attractors'][att]['phenotype'] = \
                check_cyclic_phenotype(seq)
            result['simul_result']['attractors'][att]['phenotype_data'] = seq

        elif att_type == 'point':
            att_str = state_key[att_value]
            cycstatus = check_cyclin(state_key[att_value], lbls)
            State_Apoptosis = att_str[lbls.index('State_Apoptosis')]
            if State_Apoptosis == '1':
                result['simul_result']['attractors'][att]['phenotype'] = 'A'
            elif cycstatus == '0000': 
                result['simul_result']['attractors'][att]['phenotype'] = 'Q'
            else: 
                result['simul_result']['attractors'][att]['phenotype'] = 'U'
        else: 
            result['simul_result']['attractors'][att]['phenotype'] = 'U'

    return result


def attractor_summary(data):

    A, P, Q, U = 0, 0, 0, 0

    for att in data['simul_result']['attractors']: 
        this_attr = data['simul_result']['attractors'][att]
        if this_attr['phenotype'] == 'A':
            A += this_attr['ratio'] 
        elif this_attr['phenotype'] == 'P':
            P += this_attr['ratio']
        elif this_attr['phenotype'] == 'Q':
            Q += this_attr['ratio']             
        else: 
            U += this_attr['ratio']                                     

    return A, P, Q, U

cycs = ['State_CycA','State_CycB','State_CycD','State_CycE']
