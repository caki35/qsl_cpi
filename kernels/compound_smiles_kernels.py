# -*- coding: utf-8 -*-
"""
Created on Mon May 24 00:58:48 2021

@author: Onur ÇAKI
"""

#Import Libraries
import numpy as np
import pandas as pd
import argparse
from tqdm import tqdm
import re
import os


def parse_args():
    FUNCTION_MAP = {'NLCS' : nlcs_smilarity_matrix,
                    'CLCS' : clcs_smilarity_matrix,
                    'LINGO3': lingo3_smilarity_matrix,
                    'LINGO4': lingo4_smilarity_matrix,
                    'LINGO5': lingo5_smilarity_matrix,
                    'LINGO4TF': lingo4TF_smilarity_matrix,
                    'LINGO4TFIDF': lingo4TFIDF_smilarity_matrix}
    
    
    ap = argparse.ArgumentParser()    
    ap.add_argument("smiles_path",
           help = "path of the folder in which SMILES files are located")	    
    ap.add_argument('command', choices=FUNCTION_MAP.keys())
    args = ap.parse_args()
   
    return args, FUNCTION_MAP

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def load_smile_data(data_path):    
    compound_data = {}
    for path in natural_sort(os.listdir(data_path)):
        full_path = os.path.join(data_path, path)
        if os.path.isfile(full_path) and full_path[-4:]=='.smi':
            with open(full_path) as f:
                contents = f.read()
            f.close()
            current_id = path.split('.smi')[0]
            compound_data[current_id] = contents.strip()
    return compound_data

# =============================================================================
# Similarity Calculation Functions
# =============================================================================

# Dynamic Programming implementation of NLCS problem 
def nlcs(X, Y): 
    # find the length of the strings 
    m = len(X) 
    n = len(Y) 
  
    # declaring the array for storing the dp values 
    L = [[None]*(n + 1) for i in range(m + 1)] 
  
    """Following steps build L[m + 1][n + 1] in bottom up fashion 
    Note: L[i][j] contains length of LCS of X[0..i-1] 
    and Y[0..j-1]"""
    for i in range(m + 1): 
        for j in range(n + 1): 
            if i == 0 or j == 0 : 
                L[i][j] = 0
            elif X[i-1] == Y[j-1]: 
                L[i][j] = L[i-1][j-1]+1
            else: 
                L[i][j] = max(L[i-1][j], L[i][j-1]) 
  
    # L[m][n] contains the length of LCS of X[0..n-1] & Y[0..m-1] 
    return (L[m][n]**2)/(m*n) 

#Maximal Consecutive Longest Common Subsequence starting from character 1
def MCLCS_1(SMILES1,SMILES2):    
    MCLCS1=0
    for i in range(0,min(len(SMILES1),len(SMILES2))):
        if(SMILES1[i]==SMILES2[i]):
            MCLCS1=MCLCS1+1
        else:
            break
    normolized_MCLCS1=(np.power(MCLCS1,2))/(len(SMILES1)*len(SMILES1))
    return normolized_MCLCS1

#Maximal Consecutive Longest Common Subsequence starting from character n
def MCLCS_n(X, Y): 
    
    #length of sequences
    m = len(X) 
    n = len(Y)     
    
    # LCSuff is the table with zero value initially in each cell 
   
    LCSuff = [[0 for k in range(n+1)] for l in range(m+1)] 
      
    # Initialize longest common substring 
    result = 0 
  
    # LCSuff[m+1][n+1] in bottom up fashion 
    for i in range(m + 1): 
        for j in range(n + 1): 
            if (i == 0 or j == 0): 
                LCSuff[i][j] = 0
            elif (X[i-1] == Y[j-1]): 
                LCSuff[i][j] = LCSuff[i-1][j-1] + 1
                result = max(result, LCSuff[i][j]) 
            else: 
                LCSuff[i][j] = 0
    normalized_MLCSn=(np.power(result,2))/(m*n)
    return normalized_MLCSn

#ngrams generate function
def generate_ngrams(words_list, n):
    ngrams_list = []
 
    for num in range(0, len(words_list)-(n-1)):
        ngram = ''.join(words_list[num:num + n])
        ngrams_list.append(ngram)
 
    return ngrams_list

#LINGOn funtion
def LINGOn(SMILES1,SMILES2,n):
    
    #Set ring numbers 0
    SMILES1= re.sub(r'[0-9]', r'0', SMILES1)
    SMILES2= re.sub(r'[0-9]', r'0', SMILES2)

    #generate ngrams from SMILES
    ngram1=generate_ngrams(SMILES1,n)
    ngram2=generate_ngrams(SMILES2,n)      
    ngram_merge=ngram1+ngram2      
    ngram_unique=list(set(ngram_merge)) #unique ngrams in both SMILES
      
    #Kernel Formula
    #Similarity=sum(1-|N1i-N2i/N1i+N2i|)/number of uniqiue ngram in both SMILES
    summ=0
    for lingo in ngram_unique:
        N1 = ngram1.count(lingo)
        N2 = ngram2.count(lingo)
        current_lingo = 1-(abs(N1-N2)/abs(N1+N2)) 
        summ += current_lingo
    similarity=summ/len(ngram_unique)
    return similarity

#LINGOn funtion
def TF_cos_similarity(SMILES1,SMILES2):
    
    #Set ring numbers 0
    SMILES1= re.sub(r'[0-9]', r'0', SMILES1)
    SMILES2= re.sub(r'[0-9]', r'0', SMILES2)

    #generate ngrams from SMILES
    ngram1=generate_ngrams(SMILES1,4)
    ngram2=generate_ngrams(SMILES2,4)      
    ngram_merge=ngram1+ngram2      
    ngram_unique=list(set(ngram_merge)) #unique ngrams in both SMILES
    
    #Kernel Formula
    
    #vector size
    m=len(ngram_unique)
    V_S1=np.zeros((m,1))
    V_S2=np.zeros((m,1))
    
    for i in range(0,m):
        N1 = ngram1.count(ngram_unique[i])
        N2 = ngram2.count(ngram_unique[i])
        
        if N1 == 0:
            V_S1[i] = 0
        else:
            V_S1[i] = 1 + np.log10(N1)                    
        
        if N2 == 0:
            V_S2[i] = 0
        else:
            V_S2[i] = 1 + np.log10(N2)

    similarity=np.sum(np.multiply(V_S1,V_S2))/(np.linalg.norm(V_S1)*np.linalg.norm(V_S2))
    return similarity

#LINGOn funtion
def TF_IDF_cos_similarity(SMILES1,SMILES2,dataset,n=4):
    
    #Set ring numbers 0
    SMILES1= re.sub(r'[0-9]', r'0', SMILES1)
    SMILES2= re.sub(r'[0-9]', r'0', SMILES2)
    dataset = [re.sub(r'[0-9]', r'0', SMILES) for SMILES in dataset]

    #generate ngrams from SMILES
    ngram1=generate_ngrams(SMILES1,4)
    ngram2=generate_ngrams(SMILES2,4)      
    ngram_merge=ngram1+ngram2      
    ngram_unique=list(set(ngram_merge)) #unique ngrams in both SMILES
    
    #Kernel Formula
    
    #vector size
    m=len(ngram_unique)
    V_S1=np.zeros((m,1))
    V_S2=np.zeros((m,1))
    N=len(dataset)
    
    for i in range(0,m):
        N1 = ngram1.count(ngram_unique[i])
        N2 = ngram2.count(ngram_unique[i])
        
        if N1 == 0:
            V_S1[i] = 0
        else:
            TF = 1 + np.log10(N1)
            d = [SMILES for SMILES in dataset if ngram_unique[i] in SMILES] #list of SMILES in dataset that include current LINGO
            if len(d) == 0:
                IDF = 10
            else:
                IDF=np.log10(N/len(d)) 
            V_S1[i] = TF * IDF                   
        
        if N2 == 0:
            V_S2[i] = 0
        else:
            TF = 1 + np.log10(N2)
            d = [SMILES for SMILES in dataset if ngram_unique[i] in SMILES] #list of SMILES in dataset that include current LINGO
            if len(d) == 0:
                IDF = 10
            else:
                IDF=np.log10(N/len(d))             
            V_S2[i] = TF * IDF  

    similarity=np.sum(np.multiply(V_S1,V_S2))/(np.linalg.norm(V_S1)*np.linalg.norm(V_S2))
    return similarity


# =============================================================================
# Similarity Matrix Calculation Functions
# =============================================================================
def nlcs_smilarity_matrix(compound_data):
    smiles_list = list(compound_data.values())
    id_list = list(compound_data.keys())
    similarity_matrix=np.zeros((len(smiles_list),len(smiles_list)))
    for i in tqdm(range(0,len(smiles_list))):
        for j in range(0,len(smiles_list)):
            S1=smiles_list[i]
            S2=smiles_list[j]
            similarity_matrix[i][j]=nlcs(S1,S2)
            
    sm_df = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    sm_df.to_csv('Similarity_Matrices/compound_smiles_NLCS.csv')
        
def clcs_smilarity_matrix(compound_data):
    smiles_list = list(compound_data.values())
    id_list = list(compound_data.keys())    
    similarity_matrix=np.zeros((len(smiles_list),len(smiles_list)))
    
    for i in tqdm(range(0,len(smiles_list))):
        for j in range(0,len(smiles_list)):
            S1=smiles_list[i]
            S2=smiles_list[j]
            similarity_matrix[i][j]=(MCLCS_1(S1,S2)+MCLCS_n(S1,S2)+nlcs(S1,S2))/3
            
    sm_df = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    sm_df.to_csv('Similarity_Matrices/compound_smiles_CLCS.csv')
    
def lingo4_smilarity_matrix(compound_data):
    smiles_list = list(compound_data.values())
    id_list = list(compound_data.keys()) 
    #initialize similarity matrix
    similarity_matrix=np.zeros((len(smiles_list),len(smiles_list)))
    
    #similarity calculation
    for i in tqdm(range(0,len(smiles_list))):
        for j in range(0,len(smiles_list)):
            SMILES1=smiles_list[i]
            SMILES2=smiles_list[j]
            similarity_matrix[i][j]=LINGOn(SMILES1,SMILES2,4)
            
    sm_df = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    sm_df.to_csv('Similarity_Matrices/compound_smiles_LINGO4.csv')
    
def lingo3_smilarity_matrix(compound_data):
    smiles_list = list(compound_data.values())
    id_list = list(compound_data.keys()) 
    #initialize similarity matrix
    similarity_matrix=np.zeros((len(smiles_list),len(smiles_list)))
    
    #similarity calculation
    for i in tqdm(range(0,len(smiles_list))):
        for j in range(0,len(smiles_list)):
            SMILES1=smiles_list[i]
            SMILES2=smiles_list[j]
            similarity_matrix[i][j]=LINGOn(SMILES1,SMILES2,3)
            
    sm_df = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    sm_df.to_csv('Similarity_Matrices/compound_smiles_LINGO3.csv')
    
def lingo5_smilarity_matrix(compound_data):
    smiles_list = list(compound_data.values())
    id_list = list(compound_data.keys()) 
    #initialize similarity matrix
    similarity_matrix=np.zeros((len(smiles_list),len(smiles_list)))
    
    #similarity calculation
    for i in tqdm(range(0,len(smiles_list))):
        for j in range(0,len(smiles_list)):
            SMILES1=smiles_list[i]
            SMILES2=smiles_list[j]
            similarity_matrix[i][j]=LINGOn(SMILES1,SMILES2,5)
            
    sm_df = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    sm_df.to_csv('Similarity_Matrices/compound_smiles_LINGO5.csv')

def lingo4TF_smilarity_matrix(compound_data):
    smiles_list = list(compound_data.values())
    id_list = list(compound_data.keys()) 
    
    #initialize similarity matrix
    similarity_matrix=np.zeros((len(smiles_list),len(smiles_list)))
    
    #similarity calculation
    for i in tqdm(range(0,len(smiles_list))):
        for j in range(0,len(smiles_list)):
            SMILES1=smiles_list[i]
            SMILES2=smiles_list[j]
            similarity_matrix[i][j]=TF_cos_similarity(SMILES1,SMILES2)
            
    sm_df = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    sm_df.to_csv('Similarity_Matrices/compound_smiles_LINGO4TF.csv')
    
def lingo4TFIDF_smilarity_matrix(compound_data):
    smiles_list = list(compound_data.values())
    id_list = list(compound_data.keys()) 
    
    #initialize similarity matrix
    similarity_matrix=np.zeros((len(smiles_list),len(smiles_list)))
    
    #similarity calculation
    for i in tqdm(range(0,len(smiles_list))):
        for j in range(0,len(smiles_list)):
            SMILES1=smiles_list[i]
            SMILES2=smiles_list[j]
            similarity_matrix[i][j]=TF_IDF_cos_similarity(SMILES1,SMILES2,smiles_list)
            
    sm_df = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    sm_df.to_csv('Similarity_Matrices/compound_smiles_LINGO4TFIDF.csv')    


def main():
     
    args, FUNCTION_MAP = parse_args()
        
    smiles_path = args.smiles_path
    
    compound_data = load_smile_data(smiles_path)
            
    similarity_matrix_calcation = FUNCTION_MAP[args.command]
    similarity_matrix_calcation(compound_data)

if __name__ == "__main__":
    main()