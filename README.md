# QSL_CPI
This repo provides the codes called Quasi-Supervised Learning Algorithm (QSL) for Compound-Protein Interaction Prediction that are used in following paper:
- [Quasi-Supervised Strategies for Compound-Protein Interaction Prediction](https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.202100118)

Our study based on the machine learning strategy presented in following paper:
- [Quasi-supervised learning for biomedical data analysis](https://www.sciencedirect.com/science/article/abs/pii/S0031320310002001)

If you use this repo or QSL algorithm itself in your work, please cite these two papers given above.

### Highlights of our learning algorithm:
- Compound-Protein pairs at hand are divided into two datasets (C0 and C1):
  - C1 consists of the true postive samples. (experimentally validated interactions)
  - C0 consists of unlabeled samples. (compound-protein pairs that do not have a documented interaction)
- QSL calculates posterior probability of each pair belonging to C0 and C1.
- Kolmogorov-Smirnov Method detects the threshold that draws the boundary of the overlap between the C1 and C0.
- QSL can predict the interaction profiles of unlabeled pairs accurately without requiring true negative compound-protein pairs.
- QSL is robust against data imbalance between true positives and unlabeled compound-protein pairs.
- QSL operates on similarity structure between protein and compound pairs directly without requiring a feature vector representation.

## Installation

## Download Data
We provide a script `scripts/download_data.py` for retrieving either compound or protein data from KEGG databases. You need just give a txt including id list of the data that you want to download.  
```
python scripts/download_data.py Data/Proteins/gpcr/kegg_id.txt protein
```
This script saves amino acid sequences of all given protein KEGG ids into `protein.fasta`
```
python scripts/download_data.py Data/Compounds/gpcr/kegg_id.txt compound 
```
This script saves the structure information of all given compound KEGG ids into `kegg_mol` folder in .mol file format.
