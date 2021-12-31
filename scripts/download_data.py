# -*- coding: utf-8 -*-
"""
Created on Mon May 24 00:58:48 2021

@author: Onur ÇAKI
"""

#Import Libraries
import argparse
from tqdm import tqdm
import requests # library to handle requests
import os
import urllib.request

def parse_args():
    FUNCTION_MAP = {'compound' : compound_download,
                       'protein': protein_download}
    
    ap = argparse.ArgumentParser()    
    ap.add_argument("id_list",
           help = "path of the txt files in which KEGG ids are located")	    
    ap.add_argument('command', choices=FUNCTION_MAP.keys())
    args = ap.parse_args()
   
    return args, FUNCTION_MAP

def protein_download(id_list):
    with open('protein.fasta', 'w') as file:
        for i in tqdm(range(0,len(id_list))):
            while True:
                try:
                    url="http://rest.kegg.jp/get/{}/aaseq".format(id_list[i])
                    results = requests.get(url)
                except: # You can replace Exception with something more specific.
                    continue
                else:
                    file.write("%s\n" %results.text)
                    break
    file.close()

def compound_download(id_list):
    kegg_dir = "kegg_mol"
    if not os.path.isdir("./" + kegg_dir):
        os.mkdir("./" + kegg_dir)
    for cid in tqdm(id_list):
        if not os.path.isfile("./{}/{}.mol".format(kegg_dir, cid)):
            while True:
                try:
                    url = "http://www.genome.jp/dbget-bin/www_bget?-f+m+{}".format(cid)
                    urllib.request.urlretrieve(url, "./{}/{}.mol".format(kegg_dir, cid))
                except:
                    continue
                else:
                    break
  
def main():
    args, FUNCTION_MAP = parse_args()  
    id_path = args.id_list
    with open(id_path) as f:
        id_list = f.read().splitlines()
    
    downloader = FUNCTION_MAP[args.command]
    downloader(id_list)


if __name__ == "__main__":    
    main()
