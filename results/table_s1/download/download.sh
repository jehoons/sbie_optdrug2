#/usr/bin/bash 
if [ ! -f chembl_23/chembl_23_sqlite/chembl_23.db ]; then 
    wget http://143.248.32.25/~jhsong/dataset/Compounds/ChEMBL/23/chembl_23_sqlite.tar.gz
    tar xvfz chembl_23_sqlite.tar.gz
fi

if [ ! -f chembl_23.sdf ]; then 
    wget http://143.248.32.25/~jhsong/dataset/Compounds/ChEMBL/23/chembl_23.sdf.gz
    gzip -d chembl_23.sdf.gz
fi 

if [ ! -f chembl_uniprot_mapping.txt ]; then 
    wget http://143.248.32.25/~jhsong/dataset/Compounds/ChEMBL/23/chembl_uniprot_mapping.txt
fi 

if [ ! -f CCLE_NP24.2009_Drug_data_2015.02.24.csv ]; then 
    wget http://143.248.32.25/~jhsong/dataset/Dose-Response/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv
fi

