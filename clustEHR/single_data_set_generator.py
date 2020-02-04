import pandas as pd
import os
import re
import numpy as np
import random
import data_processing.py

def generate_data(n,seed,clusters,vars, noise_ratio,  ,var_type ,seperation = None ,out_file = os.getwd() ,export = False, verbose = True ):
    """
    :param n: number of patients, ints or list of ints where the lenght is equal to number of clusters
    :param seed: int, random seed
    :param clusters: int for number of clusters or list of disease modules
    :param vars: selection of vars that a
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
               'hypertension*',
               'hypothyroidism*',
               'lung_cancer*',
               'lupus*',
               'metabolic_syndrome*',
               'osteoarthritis*',
               'osteoporosis*',
               'rheumatoid_arthritis*',
               'urinary_tract_infections*']

# sorting life out if statements
if isinstance(clusters, int) or (isinstance(clusters, list) and isinstance(clusters[0], str)):
    raise ValueError('clusters needs to be either an int representing the number of clusters or a list of diseases')
elif isinstance(clusters, int) and clusters <= len(module_list):
    clusters = random.choices(module_list, clusters)
elif isinstance(cluster, int) and clusters > len(module_list):
    raise ValueError('More clusters than possible diseasese, decrease number of clusters')
else:
    clusters = [i for i in clusters for j in module_list if i in j]
if len(clusters) == 0:
    raise ValueError('Disease names included in list were not diseases in disease module')
# generate synthea data


# read in data

# compile data

# select data

# return data