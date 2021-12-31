import pandas as pd
import requests 
import argparse
from tqdm import tqdm

def parse_args():
    ap = argparse.ArgumentParser()    
    ap.add_argument("path_of_results",
           help = "path of the csv files in which the results were saved")	    
    args = ap.parse_args()
   
    return args

def kegg_coverter(results_df):
    compound_names=[]
    protein_names=[]
    for i in tqdm(range(0,len(results_df))):
        while True:
            try:
                url_c="http://rest.kegg.jp/get/{}/".format(results_df['Compound'][i]) 
                results_c = requests.get(url_c)
                url_p="http://rest.kegg.jp/get/{}/".format(results_df['Protein'][i]) 
                results_p = requests.get(url_p)
            except: 
                continue
            else:
                current_name_c=results_c.text.splitlines()[1].replace('NAME', '').replace(';', '').strip()
                current_name_c = current_name_c.split('(')[0].strip()
                compound_names.append(current_name_c)

                current_name_p=results_p.text.splitlines()[2].replace('DEFINITION', '').replace('(RefSeq)', '').replace('NAME', '').strip()
                protein_names.append(current_name_p)
                break
            
    results_df['Protein']=protein_names
    results_df['Compound']=compound_names
    
    return results_df

def main():
    args = parse_args()  
    results_df = pd.read_csv(args.path_of_results)
    results_df = kegg_coverter(results_df)
    results_df.to_csv(args.path_of_results.split('/')[-1])  
if __name__ == "__main__":    
    main()