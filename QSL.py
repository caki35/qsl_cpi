# -*- coding: utf-8 -*-
"""  
@author: Onur Caki
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import yaml
import itertools
from tqdm import tqdm
import os
import timeit
from scripts.kegg_id_converter import kegg_coverter

# =============================================================================
# Plotting Parameters
# =============================================================================

plt.rcParams['text.color'] = '323034'
plt.rcParams['lines.markeredgecolor'] = 'black'
plt.rcParams['patch.force_edgecolor'] = True
plt.rcParams['patch.linewidth'] = 0.8
plt.rcParams['grid.color'] = 'b1afb5'

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.titlesize'] = 19
plt.rcParams['font.size'] = 16
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.framealpha'] = 0.8
plt.rcParams['legend.title_fontsize'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['legend.edgecolor'] = '0.9'
plt.rcParams['legend.borderpad'] = 0.2
plt.rcParams['legend.columnspacing'] = 1.5
plt.rcParams['legend.labelspacing'] = 0.4

plt.rcParams['axes.labelpad'] = 3
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 2


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('config', help='the alogrithm config file path')
    args = ap.parse_args()
    return args


def import_similartiy_matrix(data_path):
    if data_path.split('.')[-1] == 'txt':
        sm_matrix = pd.read_csv(data_path, index_col=0, sep='\t')
    elif data_path.split('.')[-1] == 'csv':
        sm_matrix = pd.read_csv(data_path, index_col=0, sep=',')
    else:
        print("Unvalid format. Quitting...")
        quit()

    return sm_matrix


def load_interaction_info(path):
    with open(path) as f:
        raw_info = f.read().splitlines()
    interaction_list = []
    for line in raw_info:
        protein, compound = line.split('\t')
        interaction_list.append((compound.strip(), protein.strip()))

    return interaction_list


def pair_wise_kernel(sm_compound, sm_protein, interactiondata):

    #### Modify similarity matrices so as to satisfy Mercer's Theorem ####
    from scipy import linalg

    #Drug
    Sm = sm_compound.values  # create numpy matrix
    lmbda = linalg.eigvals(Sm)  # eigen values of SM
    lmbda_min = abs(lmbda.min())  # take mimimum absolute eigen value
    identity_matrix = np.eye(len(Sm))  # create idendity matrix

    Kc = (Sm+Sm.T)/2+(identity_matrix*lmbda_min)  # Kernel Matrix for drug

    #Target
    Sm = sm_protein.values
    lmbda = linalg.eigvals(Sm)
    lmbda_min = abs(lmbda.min())
    identity_matrix = np.eye(len(Sm))

    Kp = (Sm+Sm.T)/2+(identity_matrix*lmbda_min)  # Kernel Matrix for target

    #### Combine Similarity Matrices to obtain Pair Similarities####
    #Kronocker product between two similarity matrices through numpy
    K = np.kron(Kc, Kp)
    #K = np.round(K, 5)

    compund_id = sm_compound.index.values
    protein_id = sm_protein.index.values
    merged_list = [compund_id, protein_id]
    y = []
    pairs_dict = {}
    labels_dict = {}
    i = 1
    for element in itertools.product(*merged_list):
        y.append(interactiondata.count(element))
        current_pair = {"pair{}".format(i): element}
        labels_pair = {"pair{}".format(i): interactiondata.count(element)}
        pairs_dict.update(current_pair)
        labels_dict.update(labels_pair)
        i = i+1

    pairs_pseudo = list(pairs_dict.keys())
    sm_pair = pd.DataFrame(data=K, index=pairs_pseudo, columns=pairs_pseudo)
    labels_df = pd.DataFrame(data=y, index=pairs_pseudo, columns=["label"])

    return sm_pair, labels_df, pairs_dict, labels_dict

#### Definition of functions ####


def distance_calculate(similarity_matrix, j, label_data):

    xi = similarity_matrix.iloc[:, [j]].sort_values(
        by=similarity_matrix.columns[j], ascending=False)

    labels = []
    for index in xi.index.values:
        labels.append(label_data[index])
    #append it as a new column into xi

    return labels


def posterior_calculate(di, n):
    '''
    Inputs:
    di --> the binary numpy array consist of labels of all instance 
           0. index = label of querry instance
           Other labels are sorted by similarity between the instances to which they belong and query instance
    n  --> how many isntances are used from C0 and C1 seperately to construct reference sets
    ----------------------------------------------------------------------------------------------
    Output:
    f1 --> posterior probabilty of query instance of belonging to C1 set.
    '''
    #identify k*
    k_prime = min([l[0] for l in enumerate(di) if l[1] == 1][-n],
                  [l[0] for l in enumerate(di) if l[1] == 0][-n])

    Pr_y_Ek = di[k_prime]  # P(y|k_prime-1)=1(y)
    for k in range(k_prime-1, 0, -1):
        if(di[k] == 1):  # the total number of samples with label 1 beyond this point
            l_n = di[k:].count(1)
        if(di[k] == 0):  # the total number of samples with label 0 beyond this point
            l_n = di[k:].count(0)

        Pr_y_Ek = (n/l_n)*di[k] + (1-(n/l_n))*Pr_y_Ek

    return Pr_y_Ek


def QSL_algorithm(similarity_matrix, label_data, nrange):
    '''
    it's takes three inputs:
    similarity_matrix --> Similarity values between instances
                          It must be in dataframe format. 
                          Its indices and columns name must be same
    label_data        --> It's also dataframe where the last column consist of label of instances
                          The index name, or numbers in label_data must be same with the similarity_matrix 
    nrane             --> the range of number of sample are taken each class to construct reference sets
    --------------------------------------------------------------------------
    Returns dictionary output that contains following elements:
    n_list: evaluated list of n (1,2,...,nmax)
    cost: calculated cost of which we use each n in n_list
    f0: posterior probabilities of which each samples belong to C0 when we use the n with minimum cost
    f1: posterior probabilities of which each samples belong to C1 when we use the n with minimum cost
    '''

    #number of instance
    m = len(label_data)

    #the number of instance in C1
    number_of_C1 = sum(label_data.values())

    #excract labels from dict
    instances_labels = np.fromiter(label_data.values(), dtype=int)

    #initialize posterior probabilites
    f0 = np.zeros((1, m))
    f1 = np.zeros((1, m))

    n_min = nrange[0]
    n_max = nrange[1]

    #list of the n that we will use to calculate posterior probabilities
    n_list = np.arange(n_min, n_max+1, 1)
    #initialize cost array that in which each cost for each n will keep
    cost = np.zeros((1, len(n_list)))
    for n in n_list:
        print("Calculating Posterior Probabilities for n={}".format(n))
        for i in tqdm(range(0, m)):
            # similarity between i. insantance and all other instances
            di = distance_calculate(similarity_matrix, i, label_data)
            # posterior probability of which i. instance belong C1
            f1[0][i] = posterior_calculate(di, n)
        f0 = 1-f1
        print('')
        # calculate overlap cost over the C1
        cost_C1 = 4*np.sum(np.multiply(np.multiply(f0, f1), instances_labels))
        # calculate overlap cost over the C0
        cost_C0 = 4*np.sum(np.multiply(np.multiply(f0, f1),
                           abs(instances_labels-1)))

        #normalize cost over C0 to avoid class-inbalance
        #penalize greater n to achieve generalization
        cost[0][n-1] = cost_C1+((number_of_C1/(m-number_of_C1)*cost_C0))+(2*n)

    n_list = np.reshape(n_list, (1, n_max))  # (6,) -> (1,6)

    #Find optimum n
    n_optimum = n_list[0, np.argmin(cost)]

    #by using optimum, calculate optimum posterior probabilities
    f0_optimum = np.zeros((1, m))
    f1_optimum = np.zeros((1, m))
    print("Calculating Posterior Probabilities for optimum n={}".format(n_optimum))
    for i in tqdm(range(0, m)):
        # similarity between i. insantance and all other instances
        di = distance_calculate(similarity_matrix, i, label_data)
        # posterior probabilit of which i. instance belong C1
        f1_optimum[0][i] = posterior_calculate(di, n_optimum)
    f0_optimum = 1-f1_optimum

    results = {"n_list": n_list,
               "cost": cost,
               "n_optimum": n_optimum,
               "f0": f0_optimum,
               "f1": f1_optimum}
    return results


def get_minimum(cost_dict):
    exp_dict = cost_dict.copy()
    min_n1 = min(exp_dict.keys(), key=(lambda k: exp_dict[k]))
    exp_dict.pop(min_n1, None)
    min_n2 = min(exp_dict.keys(), key=(lambda k: exp_dict[k]))
    return min_n1, min_n2


def n_list_generator(min_n1, min_n2, n_steps):
    n_list = [min_n1, min_n2]
    n_list.sort()

    while(n_list[-2]+n_steps < max(n_list)):
        n_list.append(n_list[-2]+n_steps)
        n_list.sort()

    return n_list


def QSL_algorithm_speed_up(similarity_matrix, label_data, nrange):
    '''
    Same with QSL but grid search is used to find best n with minimum cost
    '''
    import collections

    #number of instance
    m = len(label_data)

    #the number of instance in C1
    number_of_C1 = sum(label_data.values())

    #excract labels of instances names from label dataframe
    instances_labels = np.fromiter(label_data.values(), dtype=int)

    #initialize posterior probabilites
    f0 = np.zeros((1, m))
    f1 = np.zeros((1, m))

    #initialize cost array that in which each cost for each n will keep
    cost_dict = {}

    min_n1 = nrange[0]
    min_n2 = nrange[1]
    n_steps = int((min_n2-min_n1)/4)

    while (n_steps > 1):
        n_list = n_list_generator(min_n1, min_n2, n_steps)
        for n in n_list:
            if n in cost_dict.keys():
                continue
            print("Calculating Posterior Probabilities for n={}".format(n))
            for i in tqdm(range(0, m)):
                # similarity between i. insantance and all other instances
                di = distance_calculate(similarity_matrix, i, label_data)
                # posterior probability of i. instance belonging to C1
                f1[0][i] = posterior_calculate(di, n)
            print("")
            f0 = 1-f1

            # calculate overlap cost over the C1
            cost_C1 = 4 * \
                np.sum(np.multiply(np.multiply(f0, f1), instances_labels))
            # calculate overlap cost over the C0
            cost_C0 = 4 * \
                np.sum(np.multiply(np.multiply(f0, f1), abs(instances_labels-1)))

            #normalize cost over C0 to avoid class-inbalance
            #penalize greater n to achieve generalization
            cost_total = cost_C1 + \
                ((number_of_C1/(m-number_of_C1)*cost_C0))+(2*n)

            current_cost = {n: cost_total}
            cost_dict.update(current_cost)

        min_n1, min_n2 = get_minimum(cost_dict)
        n_steps = int((max(min_n2, min_n1)-min(min_n2, min_n1))/4)

    n_list = np.arange(min(min_n2, min_n1), max(min_n2, min_n1)+1, 1)
    for n in n_list:
        if n in cost_dict.keys():
            continue
        print("Calculating Posterior Probabilities for n={}".format(n))
        for i in tqdm(range(0, m)):
            # similarity between i. insantance and all other instances
            di = distance_calculate(similarity_matrix, i, label_data)
            # posterior probability of i. instance belonging to C1
            f1[0][i] = posterior_calculate(di, n)
        print("")
        f0 = 1-f1

        # calculate overlap cost over the C1
        cost_C1 = 4*np.sum(np.multiply(np.multiply(f0, f1), instances_labels))
        # calculate overlap cost over the C0
        cost_C0 = 4*np.sum(np.multiply(np.multiply(f0, f1),
                           abs(instances_labels-1)))

        #normalize cost over C0 to avoid class-inbalance
        #penalize greater n to achieve generalization
        cost_total = cost_C1+((number_of_C1/(m-number_of_C1)*cost_C0))+(2*n)

        current_cost = {n: cost_total}
        cost_dict.update(current_cost)

    #sort cost
    cost_dict = collections.OrderedDict(sorted(cost_dict.items()))

    n_list = np.fromiter(cost_dict.keys(), dtype=float)
    cost = np.fromiter(cost_dict.values(), dtype=float)

    n_list = np.reshape(n_list, (1, len(n_list)))  # (6,) -> (1,6)
    cost = np.reshape(cost, (1, len(cost)))  # (6,) -> (1,6)

    #Find optimum n
    n_optimum = min(cost_dict.keys(), key=(lambda k: cost_dict[k]))

    #by using optimum, calculate optimum posterior probabilities
    f0_optimum = np.zeros((1, m))
    f1_optimum = np.zeros((1, m))
    print("Calculating Posterior Probabilities for optimum n={}".format(n_optimum))
    for i in tqdm(range(0, m)):
        # similarity between i. insantance and all other instances
        di = distance_calculate(similarity_matrix, i, label_data)
        # posterior probabilit of which i. instance belong C1
        f1_optimum[0][i] = posterior_calculate(di, n_optimum)
    print('')
    f0_optimum = 1-f1_optimum

    results = {"n_list": n_list,
               "cost": cost,
               "n_optimum": n_optimum,
               "f0": f0_optimum,
               "f1": f1_optimum}
    return results


def QSL_algorithm_single_n(similarity_matrix, label_data, n_optimum):
    '''
    it's takes three inputs:
    similarity_matrix --> Similarity values between instances
                          It must be in dataframe format. 
                          Its indices and columns name must be same
    label_data        --> It's also dataframe where the last column consist of label of instances
                          The index name, or numbers in label_data must be same with the similarity_matrix 
    n                 --> number of sample are taken each class to construct reference sets
    --------------------------------------------------------------------------
    Returns dictionary output that contains following elements:
    n: evaluated list of n 
    f0: posterior probabilities of which each samples belong to C0 when we use the n with minimum cost
    f1: posterior probabilities of which each samples belong to C1 when we use the n with minimum cost
    '''

    #number of instance
    m = len(label_data)
    #by using optimum, calculate optimum posterior probabilities
    f0_optimum = np.zeros((1, m))
    f1_optimum = np.zeros((1, m))
    print("Calculating Posterior Probabilities for optimum n={}".format(n_optimum))
    for i in tqdm(range(0, m)):
        # similarity between i. insantance and all other instances
        di = distance_calculate(similarity_matrix, i, label_data)
        # posterior probabilit of which i. instance belong C1
        f1_optimum[0][i] = posterior_calculate(di, n_optimum)
    print('')
    f0_optimum = 1-f1_optimum

    results = {"n_optimum": n_optimum,
               "f0": f0_optimum,
               "f1": f1_optimum}
    return results


def Kolmogorov_Dmax(posteriordist_df):

    C1 = posteriordist_df[posteriordist_df["label"] == 1]
    C0 = posteriordist_df[posteriordist_df["label"] == 0]

    #copy the samples and their posterior distributions into new dataframe
    Kolmogorov_df = posteriordist_df.copy()

    #add new columns that wil be used to calculate CDF and Dmax into this dataframe
    zero_column = np.zeros((len(Kolmogorov_df), 1))
    Kolmogorov_df["F1_C0"] = zero_column  # the column of CDF for C0
    Kolmogorov_df["F1_C1"] = zero_column  # the column of CDF for C1
    # the column of difference between two CDF
    Kolmogorov_df["Dmax"] = zero_column

    #sort the posterior distributions the samples in C1 and C0 together
    sorted_df = Kolmogorov_df.sort_values(by=['f1'])

    #give equal probability of occuring to each samples in two dataset
    p0 = 1/len(C0)
    p1 = 1/len(C1)

    #inport to numpy matrix for calculations
    KSD_table = sorted_df.values

    Dmax = 0
    T = 0

    #first sample
    if KSD_table[0][0] == 0:
        KSD_table[0][2] = KSD_table[0][2]+p0
    else:
        KSD_table[0][3] = KSD_table[0][3]+p1
    KSD_table[0][4] = abs(KSD_table[0][3]-KSD_table[0][2])

    #CDF calculation
    for i in range(1, len(KSD_table)):
        if KSD_table[i][0] == 0:  # curent sample in C0
            KSD_table[i][2] = KSD_table[i-1][2]+p0
            KSD_table[i][3] = KSD_table[i-1][3]
        else:  # curent sample in C1
            KSD_table[i][3] = KSD_table[i-1][3]+p1
            KSD_table[i][2] = KSD_table[i-1][2]
        KSD_table[i][4] = abs(KSD_table[i][3]-KSD_table[i][2])

        if KSD_table[i][4] >= Dmax:  # Comparing Dmax
            Dmax = KSD_table[i][4]  # If bigger change
            F1_C0_T = KSD_table[i][2]  # Save CDF of C0
            F1_C1_T = KSD_table[i][3]  # Save CDF of C1
            T = KSD_table[i][1]  # Save current posterior as threshold

    CDF_C0 = KSD_table[:, 2]
    CDF_C1 = KSD_table[:, 3]
    posterior_samples = KSD_table[:, 1]

    ks_results = {"Dmax": Dmax,
                  "T": T,
                  "CDF_C0": CDF_C0,
                  "CDF_C1": CDF_C1,
                  "posterior_samples": posterior_samples,
                  "F1_C0_T": F1_C0_T,
                  "F1_C1_T": F1_C1_T}

    return ks_results


def plot_results(posteriordist_df, ks_results, save_dir):
    C1 = posteriordist_df[posteriordist_df["label"] == 1]
    C0 = posteriordist_df[posteriordist_df["label"] == 0]
    Dmax = ks_results['Dmax']
    T = ks_results['T']
    F1_C0 = ks_results['CDF_C0']
    F1_C1 = ks_results['CDF_C1']
    f1_all = ks_results['posterior_samples']

    F1_C0_T = ks_results['F1_C0_T']
    F1_C1_T = ks_results['F1_C1_T']

    #### Saving Results ####
    import seaborn as sns
    fig, axs = plt.subplots(2, 2)
    fig.set_figheight(8.1)
    fig.set_figwidth(10.8)

    count, bin_edges = np.histogram(posteriordist_df.iloc[:, 1])

    #create Dmax
    y_dmax = np.arange(F1_C1_T, F1_C0_T, 0.01)
    x_dmax = np.full(y_dmax.shape, T)

    sns.distplot(C1.iloc[:, 1], axlabel=False, 
                 rug=True, bins=bin_edges, ax=axs[0, 1])
    axs[0, 1].title.set_text('Distribution of f1 in C1')
    axs[0, 1].set_xlim(0, 1)
    axs[0, 1].grid(True)

    sns.distplot(C0.iloc[:, 1], axlabel=False,
                 rug=True, bins=bin_edges, ax=axs[1, 1])
    axs[1, 1].title.set_text('Distribution of f1 in C0')
    axs[1, 1].set_xlim(0, 1)
    axs[1, 1].grid(True)
    axs[1, 1].axvline(x=T, ymin=0, ymax=1, label='T', c='m')
    axs[1, 1].annotate("T= "+str(round(T, 3)), xy=(T, 5))
    axs[1, 1].set_xlabel('P(C1|x)')

    gs = axs[0, 0].get_gridspec()
    # remove the underlying axes
    for ax in axs[0:, 0]:
        ax.remove()
    axbig = fig.add_subplot(gs[0:, 0])

    axbig.plot(f1_all, F1_C0, 'r', label='C0')
    axbig.plot(f1_all, F1_C1, 'b', label='C1')
    axbig.axvline(x=T, ymin=0, ymax=1, label='T = '+str(round(T, 3)), c='m')
    axbig.plot(x_dmax, y_dmax, 'y', label='Dmax = '+str(round(Dmax, 3)))
    axbig.annotate("T = "+str(round(T, 3)), xy=(T, 0.01))
    axbig.annotate("Dmax", xy=(T, (F1_C1_T+F1_C0_T)/2))
    axbig.title.set_text('CDF of posterior probabilities (f1)')
    axbig.legend()
    axbig.set_xlim([0, 1])
    axbig.set_xlabel('P(C1|x)')
    axbig.grid(True)
    fig.savefig(os.path.join(save_dir, 'sample_distributions.png'), dpi=900)
    fig.clf()


def plot_cost_function(n_list, n_optimum, cost, save_dir):

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    plt.figure(figsize=(7, 5))
    plt.plot(n_list, cost, '-bo')
    index = np.where(n_list == n_optimum)

    plt.plot(n_optimum, cost[index[0], index[1]], '-ro')
    plt.annotate("min = {}".format(n_optimum), xy=(
        n_optimum, cost[index[0], index[1]]+2.5))
    plt.grid()
    plt.xlabel('n')
    plt.ylabel('E(n)')
    plt.title('Cost Function')
    plt.savefig(os.path.join(save_dir, 'cost_function.png'), dpi=900)


def predict_interactions(posteriordist_df, T, pairs_dict, save_dir, real_name_flag):

    estimated_prediction = posteriordist_df[(posteriordist_df["f1"] > T) & (
        posteriordist_df["label"] == 0)].sort_values(by=['f1'], ascending=False)

    drugs = []
    targets = []
    ranks = []
    for i in range(0, len(estimated_prediction)):
        (current_drug,
         current_target) = pairs_dict[estimated_prediction.index[i]]
        drugs.append(current_drug)
        targets.append(current_target)
        ranks.append("Rank {}".format(i+1))

    compounds = []
    proteins = []
    for i in range(0, len(posteriordist_df)):
        (current_compound,
         current_protein) = pairs_dict[posteriordist_df.index[i]]
        compounds.append(current_compound)
        proteins.append(current_protein)
    posteriordist_df.insert(0, "Compound", compounds, True)
    posteriordist_df.insert(1, "Protein", proteins, True)


    results = pd.DataFrame(list(zip(drugs, targets, estimated_prediction['f1'])),
                           index=ranks,
                           columns=['Compound', 'Protein', "Posterior Probability"],)

    if real_name_flag:
        print("Please wait, KEGG IDs are converted into Trivial Names")
        posteriordist_df = kegg_coverter(posteriordist_df)
        results = kegg_coverter(results)

    results.to_csv(os.path.join(save_dir, 'estimated_interaction.csv'))
    posteriordist_df.to_csv(os.path.join(
        save_dir, 'posteriors_of_all_pairs.csv'))

def main():

    args = parse_args()
    config_path = args.config
    #config_path = 'config.yml'
    with open(config_path, "r") as ymlfile:
        cfg = yaml.safe_load(ymlfile)

    #Import Similarity Imformations
    sm_compound = import_similartiy_matrix(
        cfg['qsl']['compound_similarity_matrix'])
    sm_protein = import_similartiy_matrix(
        cfg['qsl']['protein_similarity_matrix'])

    #Import Interaction Information
    interactiondata = load_interaction_info(cfg['qsl']['interaction_info'])

    sm_pair, labels_df, pairs_dict, labels_dict = pair_wise_kernel(
        sm_compound, sm_protein, interactiondata)

    if cfg['qsl']['single_n']['single_n_flag']:
        results = QSL_algorithm_single_n(sm_pair, labels_dict, cfg['qsl']['single_n']['n'])
    else:
        if cfg['qsl']['grid_search_optimizer']:
            results = QSL_algorithm_speed_up(
                sm_pair, labels_dict, cfg['qsl']['n_range'])
            print("Classification has finished. Results are preparing...")
        else:
            results = QSL_algorithm(
                sm_pair, labels_dict, cfg['qsl']['n_range'])
            print("Classification has finished. Results are preparing...")

    f1 = results['f1']
    n_optimum = results['n_optimum']
    posteriordist_df = labels_df.copy()
    posteriordist_df["f1"] = f1.T
    if not os.path.exists(cfg['qsl']['save_dir']):
        os.mkdir(cfg['qsl']['save_dir'])

    #Read csv
    #posteriordist_df = pd.read_csv('NR_MCS/posteriors_of_all_pairs.csv', index_col = 0)

    #### Kolmogorov-Smirnov Test ####
    ks_results = Kolmogorov_Dmax(posteriordist_df)

    #export outputs
    plot_results(posteriordist_df, ks_results, cfg['qsl']['save_dir'])
    if not cfg['qsl']['single_n']['single_n_flag']:
        cost = results['cost']
        n_list = results['n_list']
        plot_cost_function(n_list, n_optimum, cost, cfg['qsl']['save_dir'])
    predict_interactions(
        posteriordist_df, ks_results['T'], pairs_dict, cfg['qsl']['save_dir'], cfg['qsl']['trivial_names'])

    print("Finished")

if __name__ == "__main__":
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    total_time = stop - start

    # output running time in a nice format.
    mins, secs = divmod(total_time, 60)
    hours, mins = divmod(mins, 60)
    print("")
    print("Total running time: %d:%d:%d.\n" % (hours, mins, secs))
