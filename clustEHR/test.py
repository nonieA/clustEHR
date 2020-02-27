import single_data_set_generator as sg
import numpy as np
import pandas as pd
import re
import os
import noise_seperation as ns
import csv
import data_generation as dg
import seaborn as sb
from statistics import stdev

file_list = os.listdir('C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\')

df = sg.generate_data(n=1000, seed=11, clusters=21, vars=None,
                      noise_var_ratio=[[1, 1, 1], [0, 0, 0]], var_n=20,
                      out_file='C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\')

file_list2 = os.listdir('C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\')
file = [i for i in file_list2 if i not in file_list and 'output' in i]

df_y = pd.read_csv('C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\' + file + '\\labels.csv')
var_selector = pd.read_csv('C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\' + file + '\\varstypes.csv')
df = pd.merge(df, df_y,how = 'inner', on = 'PATIENT')
df = df.drop('Unnamed: 0', axis = 1)



def seps(df,disease1,disease2,col,vars,perc):
    s1 = df.loc[df['DISEASE'] == disease1, col]
    s2 = df.loc[df['DISEASE'] == disease2, col]
    if vars.at[vars.index[vars['vars'] == col][0],'type'] == 'cont':
        return(ns.percentile_dif(s1,s2,perc))
    else:
        return(ns.mode_dif(s1,s2))

col_list = var_selector.loc[var_selector['noise'] == 'feature','vars'].tolist()
disease_list = df['DISEASE'].unique()


def mean_finder(perc):
    test_dict = {i + '_' + j: [seps(df, i, j, k, var_selector, perc) for k in col_list] for i in disease_list for j in
                 disease_list if i != j}

    test_df = pd.DataFrame(test_dict)
    return(test_df)

perc_dict = {str(i): mean_finder(i) for i in [5,10,20,25,40]}

def sep_getter(s1):
    return(s1.dropna().sum()/len(s1.dropna()))

total_dict = {k: list(map(lambda x: sep_getter(df[x]),df.columns)) for k,df in perc_dict.items()}
total_df = pd.DataFrame(total_dict)

for i in total_df.columns:
    sb.distplot(total_df[i], hist = False, label = i)

total_df.index = perc_dict['5'].columns

std_list = [stdev(perc_dict['10'][i].dropna()) for i in perc_dict['10'].columns]

cut_offs = pd.DataFrame({'Disease_comb':perc_dict['5'].columns, 'Means': total_df['10'],'stddev':std_list})
cut_offs = (cut_offs.assign(Sep_group= pd.qcut(cut_offs['Means'],5,labels = False),
                            Hom_group= pd.qcut(cut_offs['stddev'], 3, labels = False)))

cut_offs['disease1'] = [re.sub('_.*', '', i) for i in cut_offs['Disease_comb'].tolist()]
cut_offs['disease2'] = [re.sub('.*_', '', i) for i in cut_offs['Disease_comb'].tolist()]
sep_df = cut_offs.pivot(index = 'disease1', columns = 'disease2', values = 'Sep_group').rename_axis(None)
sep_df = sep_df[sep_df.index.tolist()]
sep_mat = sep_df.to_numpy()
std_df = cut_offs.pivot(index = 'disease1', columns = 'disease2', values = 'Hom_group').rename_axis(None)
std_df = std_df[std_df.index.tolist()]

order = random.sample(range(len(sep_mat)),len(sep_mat))

icd_10 = pd.read_csv('C:\\Users\\nonie\\Documents\\icd10cm_codes_2019.txt', sep = '\t')
icd_list = icd_10.iloc[:,0].to_list()
module_list = [re.sub('\*','',i) for i in module_list]
module_list = [re.sub('_',' ',i) for i in module_list]
module_list[module_list == 'urinary tract infections'] = 'urinary tract infection'
def disease_checker(str1,str2):
     str_check = str1.lower()
     list_in = str2.split()
     if all(i in str_check for i in list_in):
         return(True)
     else:
         return(False)


codes = [[i,j] for i in icd_list for j in module_list if disease_checker(i,j)]
icd_df = pd.DataFrame({'disease':[i[1] for i in codes],
                       'codes':[re.sub(' .*','',i[0]) for i in codes],
                       'description':[re.sub('\A\S* *','',i[0]) for i in codes]})
icd_df = icd_df.sort_values('disease')
icd_df = icd_df[~icd_df['codes'].str.contains('Z')]
icd_df = icd_df[icd_df['codes'].str.len() <= 4]
icd_df['codes'] = icd_df['codes'].str.slice(0,2)
icd_df = icd_df.drop_duplicates(['disease','codes'])
#Emergencany get data

file_list = os.listdir('C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\')

def csv_read_count(file):
    if os.path.isfile('C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\' + file + '\\patstest.txt'):
        row_count = sum(1 for row in csv.reader(open('C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\' + file + '\\patstest.txt')))
    else:
        row_count = 0
    if row_count >= 1000:
        return(file)

new_file_list = list(map(csv_read_count,file_list))
new_file_list = [i for i in new_file_list if i]
new_file_list.remove('gallstones_2020-02-22_11')
get_df_list = [i for i in new_file_list
               if os.path.isdir('C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\'+ i + '\\csv')]

imp_folder = 'C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\'

def get_diseases(file_name):
    out_folder = 'C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\'
    file_list = dg._read_files(file_name, n=1000, out_file=out_folder)
    disease_cut = re.sub('_2.*', '',file_name)
    df = dg.full_out(disease_cut, df_list=file_list, write_out=(out_folder + file_name + '/'))
    dg._remove_files(out_folder + file_name + '/')
    return (df)

df_list = [get_diseases(i) for i in get_df_list]

def read_thing(file):
    ds = re.sub('_2.*', '', file)
    if os.path.isfile('C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\' + file + '\\' + ds + 'clean.csv'):
        df = pd.read_csv('C:\\Users\\nonie\\Documents\\clustEHR\\run_through\\' + file + '\\' + ds + 'clean.csv')
        return(df)

disease_list = [read_thing(file) for file in new_file_list]

disease_list = [df.drop('Unnamed: 0', axis = 1) for df in disease_list]

disease_list = disease_list + df_list