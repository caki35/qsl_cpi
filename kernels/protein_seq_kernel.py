# -*- coding: utf-8 -*-
"""
Created on Thu May 27 21:11:48 2021

@author: Onur Caki
"""

#Import fundamental Libraries
import numpy as np
import pandas as pd
import time
from Bio import SeqIO
import argparse
from tqdm import tqdm
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import MutableSeq


#Create simith-waterman aligner with default values of EMBOSS Water
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.substitution_matrix = substitution_matrices.load("blosum62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5


def parse_args():
    ap = argparse.ArgumentParser()    
    ap.add_argument("data_path",
           help = "path of the .fast file in which the amino acid sequences are located")	
    args = ap.parse_args()
   
    return args

def normalized_smith_waterman(Seq1,Seq2):
    
    ## Change amino aicd U with unkown X in order to perform simith-waterman alignment
    mutable_seq=MutableSeq(str(Seq1))
    if 'U' in mutable_seq:
        for i in range(0,len(mutable_seq)):
            if mutable_seq[i]=="U":
                mutable_seq[i] ="X"
        Seq1.seq= mutable_seq.toseq()
    
    mutable_seq=MutableSeq(str(Seq2))
    if 'U' in mutable_seq:
        for i in range(0,len(mutable_seq)):
            if mutable_seq[i]=="U":
                mutable_seq[i] ="X"
        Seq2.seq = mutable_seq.toseq()
    
    ## Calculate normalized simith-waterman score
    SW_ij = aligner.score(Seq1,Seq2)
    SW_ii = aligner.score(Seq1,Seq1)
    SW_jj = aligner.score(Seq2,Seq2)
    similarity_score = SW_ij/np.sqrt(SW_ii*SW_jj)   
    return similarity_score

def similarity_matrix_calculation(aminoacid_data):
    similarity_matrix=np.zeros((len(aminoacid_data),len(aminoacid_data)))
    for i in tqdm(range(0,len(aminoacid_data))):
        for j in tqdm(range(0,len(aminoacid_data)),leave=False):
            similarity_matrix[i][j]=normalized_smith_waterman(aminoacid_data[i].seq, aminoacid_data[j].seq)
            time.sleep(0.01)
    return similarity_matrix

def main():
    args = parse_args()  
    data_path = args.data_path
    with open(data_path) as file:
        aminoacid_data = list(SeqIO.parse(file,'fasta'))
    file.close()

    id_list=[]
    for i in range(0,len(aminoacid_data)):
        id_list.append(aminoacid_data[i].id)
        
    similarity_matrix = similarity_matrix_calculation(aminoacid_data)
    
    sm_mcs_rdkit = pd.DataFrame(data=similarity_matrix, index=id_list, columns=id_list)
    
    sm_mcs_rdkit.to_csv('Similarity_Matrices/protein_NSWA.csv')


if __name__ == "__main__":
    main()



