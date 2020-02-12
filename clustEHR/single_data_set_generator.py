import pandas as pd
import os
import re
import numpy as np
import random
import data_generation as dg
import data_processing as dp
import noise_seperation as ns
import datetime as dt
import warnings
def generate_data(n, seed, clusters, vars, noise_var_ratio, var_n,
                  seperation = None, priority = var_n , out_file = os.getwd(), export = False, verbose = True ):
    """
    :param n: number of patients, ints or list of ints where the lenght is equal to number of clusters
    :param seed: int, random seed
    :param clusters: int for number of clusters or list of disease modules
    :param vars: selection of vars to return
    :param noise_ratio: ratio of noise variables to every 1 feature
    :param seperation: seperation of clusters as defined by gaussian distribution
    :param out_file: file to export data (default: working directory)
    :param export: export data at all points (default: false)
    :param verbose: return updates of funcitons
    :return: 1 cluster test data set
    """
    module_list = ['appendicitis*',
                   'asthma*',
                   'atopy*',
                   'breast_cancer*',
                   'bronchitis*',
                   'colorectal_cancer*',
                   'copd',
                   'cystic_fibrosis*',
                   'dementia',
                   'dermatitis*',
                   'epilepsy*',
                   'fibromyalgia*',
                   'food_allergies*',
                   'gallstones*',
                   'gout*',
                   'hypothyroidism*',
                   'lung_cancer*',
                   'lupus*',
                   'metabolic_syndrome*',
                   'osteoarthritis*',
                   'osteoporosis*',
                   'rheumatoid_arthritis*',
                   'urinary_tract_infections*']

    # sorting life out if statements
    if not isinstance(clusters, int) or (isinstance(clusters, list) and isinstance(clusters[0], str)):
        raise ValueError('clusters needs to be either an int representing the number of clusters or a list of diseases')
    elif isinstance(clusters, int) and clusters <= len(module_list):
        clusters = random.sample(module_list, clusters)
    elif isinstance(cluster, int) and clusters > len(module_list):
        raise ValueError('More clusters than possible diseasese, decrease number of clusters')
    else:
        clusters = [i for i in clusters for j in module_list if i in j]

    if len(clusters) == 0:
        raise ValueError('Disease names included in list were not diseases in disease module')

    if isinstance(n,int):
         n = [n for i in range(len(clusters))]

    if noise_var_ratio == None:
        noise_var_ratio = [[3,3],[1,1]]

    # generate synthea data

    def one_dis(n, disease, seed, out_folder):
        dg._disease_counter(n,disease, seed, out_folder)
        file_name = (re.sub('\*',"",disease, count = 0) +
                    "_" +
                    str(dt.datetime.now().date()) +
                    "_" +
                    str(seed)
                    )

        file_list = dg._read_files(file_name, out_file= out_folder)
        disease_cut = re.sub('\*','',disease)
        df = dg.full_out(disease_cut, df_list = file_list, write_out = (out_folder + file_name + '/'))
        dg._remove_files(out_file + file_name + '/')
        return(df)

    disease_list = [one_dis(n[i], clusters[i], seed, out_folder = out_file) for i in range(len(n))]


    comb_df, excepts = dp._combine_disease_dfs(disease_list)
    df_X, df_y, outcomes = dp.data_clean(comb_df, comb = excepts)

    if vars == None:
        var_selector = pd.merge(ns.rf_noise(df_X, df_y),ns.var_type(df_X), how = 'left', on = 'vars')
        var_count = dp.var_ratio_returner(var_selector,var_n = var_n, noise_var_ratio=noise_var_ratio, priority = priority)
        vars = dp.var_getter(var_count,var_selector)

    df_fin = df_X[vars]


























    # read in data

    # compile data

    # select data

    # return data