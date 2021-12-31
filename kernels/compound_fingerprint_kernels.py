# -*- coding: utf-8 -*-
"""
Created on Sun May 30 00:15:58 2021

@author: Onur Çakı
"""

import numpy as np
import pandas as pd
import argparse
import os
from tqdm import tqdm
import re
from rdkit import Chem
from rdkit.Chem import AllChem #Library for Morgan/Circular Fingerprint
from rdkit.Chem import MACCSkeys
try:
    import kcfconvoy as kcf
except Exception as e:
    print("WARNING", e)
    print("You can't use KCFS")

def parse_args():
    FINGERPRINT_MAP = {'KCFS' : kcfs_similarity_matrix,
                       'ECFP4': ECFP4_similarity_matrix,
                       'FCFP4': FCFP4_similarity_matrix,
                       'MACCS': MACCS_similarity_matrix}
    
    ap = argparse.ArgumentParser()    
    ap.add_argument("data_path",
           help = "path of the .fast file in which the amino acid sequences are located")
    ap.add_argument('command', choices=FINGERPRINT_MAP.keys())	
    args = ap.parse_args()
   
    return args, FINGERPRINT_MAP

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def load_data(path):
    mol_files = []
    id_list = []
    for maindir, subdir, file_name_list in os.walk(path):
        for filename in file_name_list:
            apath = os.path.join(maindir, filename)
            if apath[-4:]=='.mol':
                mol_files.append(apath)
                id_list.append(filename.split('.mol')[0])
    return natural_sort(mol_files), natural_sort(id_list)


# =============================================================================
# Similarity Calculation Functions
# =============================================================================
def kcfs_similarity(kcf_vec_1, kcf_vec_2):
    kegg_atom_levels = set(["atom_species", "atom_class", "kegg_atom"])
    l_count_1 = []
    l_count_2 = []
    n_nodes=list(range(99))
    for ele, label in kcf_vec_1.kcf_vec.items():
        if not label["n_nodes"] in n_nodes:
            continue
        if not label["ele_level"] in kegg_atom_levels:
            continue
        if ele in kcf_vec_2.kcf_vec.keys():
            l_count_1.append(kcf_vec_1.kcf_vec[ele]["count"])
            l_count_2.append(kcf_vec_2.kcf_vec[ele]["count"])
        else:
            l_count_1.append(kcf_vec_1.kcf_vec[ele]["count"])
            l_count_2.append(0)
    for ele, label in kcf_vec_2.kcf_vec.items():
        if not label["n_nodes"] in n_nodes:
            continue
        if not label["ele_level"] in kegg_atom_levels:
            continue
        if ele not in kcf_vec_1.kcf_vec.keys():
            l_count_1.append(0)
            l_count_2.append(kcf_vec_2.kcf_vec[ele]["count"])

    np_count_1 = np.array(l_count_1)
    np_count_2 = np.array(l_count_2)

    only_1 = np.sum(np.fmax(0, (np_count_1 - np_count_2)))
    only_2 = np.sum(np.fmax(0, (np_count_2 - np_count_1)))
    both_12 = np.sum((np.minimum(np_count_1, np_count_2)))

    if only_1 + only_2 + both_12 == 0:
        x = 0
    else:
        x = both_12 / (only_1 + only_2 + both_12)

    return x 

def ecfp4_similarity(mol1, mol2):
    fps1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024) 
    fps2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024) 
    V1=np.array(fps1)
    V2=np.array(fps2)
    similarity=float(np.dot(V1,V2)/(np.linalg.norm(V1)*np.linalg.norm(V2))) #cosine Similarity
    return similarity

def fcfp4_similarity(mol1, mol2):
    fps1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024, useFeatures=True)  
    fps2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024, useFeatures=True) 
    V1=np.array(fps1)
    V2=np.array(fps2)
    similarity=float(np.dot(V1,V2)/(np.linalg.norm(V1)*np.linalg.norm(V2))) #cosine Similarity
    return similarity
    
def maccs_similarity(mol1, mol2):
    fps1 = MACCSkeys.GenMACCSKeys(mol1)  
    fps2 = MACCSkeys.GenMACCSKeys(mol2) 
    V1=np.array(fps1)
    V2=np.array(fps2)
    similarity=float(np.dot(V1,V2)/(np.linalg.norm(V1)*np.linalg.norm(V2))) #cosine Similarity
    return similarity
    
# =============================================================================
# Similarity Matrix Calculation Functions
# =============================================================================

# KCF-S --> KEGG Chemical Function and Substructures. 
def kcfs_similarity_matrix(mol_files, id_list):
    kcfmat = kcf.KCFmat()
    for i in tqdm(range(0,len(mol_files))):
        kcfmat.input_molfile(mol_files[i],id_list[i])
        
    similarity_matrix=np.zeros((len(mol_files),len(mol_files)))

    for i in tqdm(range(0,similarity_matrix.shape[0])):
        for j in range(0,similarity_matrix.shape[1]):
            similarity_matrix[i][j]=kcfs_similarity(kcfmat.kcf_vecs[i], kcfmat.kcf_vecs[j])
            
    sm_df = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    sm_df.to_csv('Similarity_Matrices/compound_fp_KCFS.csv')  
    
# ECFP --> Extended-Connectivity Fingerprints
def ECFP4_similarity_matrix(mol_files, id_list):
    suppl=[]
    for mol in mol_files:
        m=Chem.MolFromMolFile(mol)
        suppl.append(m) 
        
    similarity_matrix=np.zeros((len(mol_files),len(mol_files)))

    for i in tqdm(range(0,similarity_matrix.shape[0])):
        for j in range(0,similarity_matrix.shape[1]):
            similarity_matrix[i][j]=ecfp4_similarity(suppl[i], suppl[j])
            
    sm_df = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    sm_df.to_csv('Similarity_Matrices/compound_fp_ECFP4.csv')  

# FCFP --> Functional-Class Fingerprints. 
def FCFP4_similarity_matrix(mol_files, id_list):
    suppl=[]
    for mol in mol_files:
        m=Chem.MolFromMolFile(mol)
        suppl.append(m) 
        
    similarity_matrix=np.zeros((len(mol_files),len(mol_files)))

    for i in tqdm(range(0,similarity_matrix.shape[0])):
        for j in range(0,similarity_matrix.shape[1]):
            similarity_matrix[i][j]=fcfp4_similarity(suppl[i], suppl[j])
            
    sm_df = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    sm_df.to_csv('Similarity_Matrices/compound_fp_FCFP4.csv') 
    
# The MACCS (Molecular ACCess System)
def MACCS_similarity_matrix(mol_files, id_list):
    suppl=[]
    for mol in mol_files:
        m=Chem.MolFromMolFile(mol)
        suppl.append(m) 
        
    similarity_matrix=np.zeros((len(mol_files),len(mol_files)))

    for i in tqdm(range(0,similarity_matrix.shape[0])):
        for j in range(0,similarity_matrix.shape[1]):
            similarity_matrix[i][j]=maccs_similarity(suppl[i], suppl[j])
            
    sm_df = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    sm_df.to_csv('Similarity_Matrices/compound_fp_MACCS.csv')  
    

def main():
    args, FINGERPRINT_MAP = parse_args()  
    data_path = args.data_path
    mol_files, id_list = load_data(data_path)
    
    similarity_matrix_calcation = FINGERPRINT_MAP[args.command]
    similarity_matrix_calcation(mol_files, id_list)


if __name__ == "__main__":
    main()