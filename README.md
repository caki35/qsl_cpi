# QSL_CPI
This repo provides the codes used in the study on Quasi-Supervised Learning Algorithm (QSL) for Compound-Protein Interaction Prediction. Please check the following link for a more detailed explanation of this work.
- [Quasi-Supervised Strategies for Compound-Protein Interaction Prediction](https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.202100118)

Our study is based on the machine learning strategy presented in the following paper:
- [Quasi-supervised learning for biomedical data analysis](https://www.sciencedirect.com/science/article/abs/pii/S0031320310002001)

If you use this repo or QSL algorithm itself in your work, please cite these two papers given above.

### Highlights of our method:
- Compound-Protein pairs at hand are divided into two datasets:
  - C1 consists of the true postive samples. (experimentally validated interactions)
  - C0 consists of unlabeled samples. (compound-protein pairs that do not have a documented interaction)
- QSL calculates posterior probability of each pair belonging to C0 and C1.
- Kolmogorov-Smirnov Method detects the threshold that draws the boundary of the overlap between the C1 and C0.
- QSL can predict the interaction profiles of unlabeled pairs accurately without requiring true negative compound-protein pairs.
- QSL is robust against data imbalance between true positives and unlabeled compound-protein pairs.
- QSL operates on similarity structure between protein and compound pairs directly without requiring a feature vector representation.

## Installation
'''bash
git clone https://github.com/caki35/qsl_cpi.git
cd qsl_cpi
'''
Istallation using conda 
```bash
conda create --name qsl_cpi --file requirements_conda_env.txt
conda activate qsl_cpi
```
Istallation using pip 
```bash
pip install requirements.txt
```

## Download Data
We provide a script `scripts/download_data.py` for retrieving either compound or protein data from KEGG databases. You need just give a txt file including id list of the data that you want to download.  
```bash
python scripts/download_data.py Data/Proteins/gpcr/kegg_id.txt protein
```
This script saves amino acid sequences of all given protein KEGG ids into a unified file named `protein.fasta`.
```bash
python scripts/download_data.py Data/Compounds/gpcr/kegg_id.txt compound 
```
This script saves the structure information of all given compound KEGG ids into `kegg_mol` folder in .mol file format.

### Generating SMILES strings
You can convert mol files of compounds into SMILES representations using `scripts/smiles_converter.sh`  
It requires [molconverter console program of JCHEM](https://chemaxon.com/products/jchem-engines/download#jchem). The program defines SMILES by following
[Daylight’s SMILES specification rules](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html). You must install this program and provide its sysytem path into `scripts/smiles_converter.sh`
```bash
bash scripts/smiles_converter.sh MOL_FOLDER_PATH
```
All .mol files in given folder are converted into SMILES representation and saved into folder `smiles_output`.

## Similarity Calculation 
To use QSL for CPI prediction, you must construct similarity matrices for both compound and protein data first. We provide various options to do so. If you have your own similarity matrices, you can pass this step.
### Compound-Compound Similarity Measure

**SMILES Kernels**  
To calculate similarity between SMILES pairs, run following script with folder path including SMILES files. You must choose one of the kernel options:
- 'NLCS' : Normalized Longest Common Subsequence (NLCS)
- 'CLCS' : Combination of Longest Common Subsequence Models (CLCS)
- 'LINGO3': LINGO-3 Similarity
- 'LINGO4': LINGO-4 Similarity
- 'LINGO5': LINGO-5 Similarity
- 'LINGO4TF': LINGO-4 Based Term Frequency (TF) Cosine Similarity
- 'LINGO4TFIDF': Frequency-Inverse Document Frequency (TF-IDF) Cosine Similarity

```bash
python kernels/compound_smiles_kernels.py SMILES_FOLDER_PATH LINGO4TF
```
**Graph Kernel**  
To calculate similarity between compound graphs in .mol file, run following script with folder path including .mol files.
```bash
python kernels/compound_graph_kernel.py Data/Compounds/GPCR/GPCR_kegg_mol
```
**Molecular-Fingerprints Kernel**  
To calculate similarity between compounds using their fingerprints, run following script with folder path including .mol files. You must choose one of the fingerprint options:
- 'KCFS' : KEGG Chemical Function and Substructures
- 'ECFP4': Extended-Connectivity Fingerprints
- 'FCFP4': FunctionalClass Fingerprintsx
- 'MACCS': Molecular ACCess System

```bash
python kernels/compound_fingerprint_kernels.py MOL_FOLDER_PATH ECFP4
```

The results will be saved into `Similarity_Matrices` in .csv file format.

### Protein-Protein Similarity Measure
Run the following script to construct a similarity matrix for proteins using Smith-Waterman Algorithm. You need to input`.fasta` file that contains amino acid sequences of all proteins.
```bash
python kernels/protein_seq_kernel.py Data/Proteins/GPCR/GPCR_protein.fasta
```
## Run QSL algorithm
```bash
python QSL.py CONFIG_FILE_PATH
```
The detailed of `config.yml`:  
`save_dir`:The directory in which to save the output of the QSL algorithm.  
`compound_similarity_matrix`: The path of similarity matrix for compounds.  
`protein_similarity_matrix`: The path of similarity matrix for proteins.  
`interaction_info`: The path of .txt file including interaction information at hand.  
`n_range`: The range of *n* which is the main parameter that QSL optimizes.  
`grid_search_optimizer`: if **True**, the algorithm performs grid search during finding the best *n* instead of considering each *n* in range.  
`single_n_flag`: if **True**, QSL algorithm runs at a certain  *n* instead of optimizing it. 

```yaml
qsl:
  save_dir: results_TFIDF_NR #the directory to which the results are saved
  compound_similarity_matrix: Similarity_Matrices/Nuclear_Receptor/compound_smiles_LINGO4TFIDF.csv
  protein_similarity_matrix: Similarity_Matrices/Nuclear_Receptor/protein_NSWA_nr.csv
  interaction_info: Data/interaction_data_nr.txt
  n_range: [1,50]
  grid_search_optimizer: False
  single_n:
    single_n_flag: False
    n: 21
```

The result folder will contain 4 files:
- The graph of cost function with which **n** is optimized.
- Posterior probability distributions of compound-protein pairs to belong to the true positive dataset and their Kolmogorov-Smirnov Analysis.
- The list of posterior probabilities of all pairs.
- The list of predicted interactions.
![image](https://user-images.githubusercontent.com/61325251/147421501-055b9700-4f97-4e58-b59c-931fea9cc98b.png)

![sample_distributions](https://user-images.githubusercontent.com/61325251/147421589-4ed34eb2-7a82-40db-9e3f-ccc98378d5c1.png)

Cost function when grid search open.
![cost_function](https://user-images.githubusercontent.com/61325251/147421598-577f7fd6-3290-467f-b54b-02aac0c9eefb.png)

Cost function when grid search close.
![cost_function](https://user-images.githubusercontent.com/61325251/147422122-70dfe22d-f464-41eb-a194-e473cb6c6f0d.png)



