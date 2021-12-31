# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 04:57:43 2021

@author: Onur Çakı
"""

import numpy as np
import pandas as pd
import re
from rdkit import Chem
from rdkit.Chem import rdFMCS
import argparse
import os
from tqdm import tqdm
import time

def parse_args():
    ap = argparse.ArgumentParser()    
    ap.add_argument("data_path",
           help = "path of the folder in which mol files are located")	
    ap.add_argument("-ds", "--dice_similarity", action='store_true',
       help = "Use dice similarity")
    args = ap.parse_args()
   
    return args
	
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def load_data(data_path):    
    mol_files = []
    id_list = []
    for path in os.listdir(data_path):
        full_path = os.path.join(data_path, path)
        if os.path.isfile(full_path) and full_path[-4:]=='.mol':
            mol_files.append(Chem.MolFromMolFile(full_path))
            id_list.append(path.split('.mol')[0])
            
    return mol_files, id_list
	
def jaccard_similarity(mol1,mol2):
    
    mols=[mol1,mol2] 
    res=rdFMCS.FindMCS(mols, bondCompare=rdFMCS.BondCompare.CompareOrderExact, matchValences=True, ringMatchesRingOnly=True,completeRingsOnly=True) #find MCS
    G12=float(res.numAtoms+res.numBonds)         #cardinality of MCS
    G1=float(mol1.GetNumAtoms()+mol1.GetNumBonds())    #cardinality of G1
    G2=float(mol2.GetNumAtoms()+mol2.GetNumBonds())    #cardinality of G2
    similarity_score = float(G12/(G1+G2-G12))
    
    return similarity_score

def dice_similarity(mol1,mol2):
    
    mols=[mol1,mol2] 
    res=rdFMCS.FindMCS(mols, bondCompare=rdFMCS.BondCompare.CompareOrderExact, matchValences=True, ringMatchesRingOnly=True,completeRingsOnly=True) #find MCS
    G12=float(res.numAtoms+res.numBonds)         #cardinality of MCS
    G1=float(mol1.GetNumAtoms()+mol1.GetNumBonds())    #cardinality of G1
    G2=float(mol2.GetNumAtoms()+mol2.GetNumBonds())    #cardinality of G2    
    similarity_score = float(G12**2/(G1*G2))
    
    return similarity_score

def similarity_calculation(mol_files,flag):
    #initialize similarity matrix
    similarity_matrix=np.zeros((len(mol_files),len(mol_files)))
    
    for m1 in tqdm(range(0,len(mol_files))):
        for m2 in tqdm(range(0,len(mol_files)),leave=False):
            mol1 = mol_files[m1] #take i. mol
            mol2 = mol_files[m2] #take j. mol
            if flag:
                similarity_score = dice_similarity(mol1,mol2)
            else:
                similarity_score = jaccard_similarity(mol1,mol2)
                similarity_matrix[m1][m2]= similarity_score
            time.sleep(0.01)

            
    return similarity_matrix
                    
def main():
    args = parse_args()  
    
    mol_files, id_list = load_data(args.data_path)
    
    similarity_matrix = similarity_calculation(mol_files,args.dice_similarity)
    
    sm_mcs_rdkit = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    
    sm_mcs_rdkit.to_csv('Similarity_Matrices/compound_MCSjaccard.csv')

if __name__ == "__main__":
    main()


