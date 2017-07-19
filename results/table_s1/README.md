### Table S1. Compound similarity analysis

#### (**A**) Prepare dataset from ChEMBL

*1. Filt subset compounds with target information*

*2. Extract uniprot names*

*3. Convert uniprot names into gene names by using uniprot.org*

*4. Prepare compound list*

#### (**B**) Calculate similarity 

*1. Collect drug information*

*2. Compute drug-drug similarity matrix*

*3. Labels for the similarity matrix*

*4. Prepare compound-target dictionary*

#### (**C**) CCLE compound analysis 

여기서 CCLE 데이터(약물, 도즈반응커브)를 활용함으로써 시뮬레이션 데이터를 검증하고자 한다. 그러면, explicit target information 과 implicit target information을 이용하는 것이 가능하다. 

*1. Prepare treatment data (CCLE)*

First, extract ccle drugs from ccle dataset, and then curate chembl_id from chembl database site. 

*2. Infer alternative targets based on compound similarity*

We are using CCLE to validate network models. Thus, we need to infer alternative targets of CCLE compounds. 