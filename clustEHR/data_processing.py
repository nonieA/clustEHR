"""
Takes multiple disease DFs and processes them to make one clusterable dataset
"""
import json
import re

import pandas as pd
import os
import numpy as np
from sklearn.impute import KNNImputer
from sklearn.preprocessing import OneHotEncoder
import random
import warnings

def _combine_disease_dfs(disease_list,description=False):
    """
    takes individiul diseases and makes them one dataset
    :param disease_list: list of pd.dfs or file paths
    :return: one df and concat dict
    """
    def list_checker(x):
        if isinstance(x, pd.DataFrame):
            return(x)
        elif os.path.isfile(x) and '.csv' in  x:
            return(pd.read_csv(x).drop('Unnamed: 0', axis = 1))
        else:
            raise ValueError('There is a problem in that the input is neither a df or a directory to a df so please can you change that')

    df_disease_list = [list_checker(i) for i in disease_list]

    if description:
        df = df_disease_list[0]
    else:
        df = pd.concat(df_disease_list, sort=False)
        for i in df['DISEASE'].unique().tolist():
            if i in df.columns:
                df.drop(columns = i)
    df = df.replace([np.inf, -np.inf], np.nan)

    return df

def find_cols(cat,var,col_dict_list,df):
    col_list = [set(i[cat]) for i in col_dict_list]
    col_list = [j for i in col_list for j in i]
    df_break = [df[df['DISEASE'] == i] for i in df['DISEASE'].unique()]
    if ((col_list.count(var) >= len(col_dict_list)) and all(len(i[var].dropna()) > 0 for i in df_break)):
        return True
    else:
        return False


def imputer(imp_df):
    imp = KNNImputer()
    out_df = imp.fit_transform(imp_df)
    out_df = pd.DataFrame(out_df, columns=imp_df.columns)
    return out_df

def get_cat_df(cat,description,col_dict_list,df):
    if description:
        cat_cols = [j for i in col_dict_list for j in i[cat] if (len(df[df[j] > 0]) / len(df)) >= 0.4]
    else:
        cat_cols = set([j for i in col_dict_list for j in i[cat] if (find_cols(cat,j,col_dict_list,df)and((len(df[df[j]>0])/len(df))>= 0.4))])
    cat_df = df[cat_cols]
    cat_df = cat_df.fillna(0)
    return cat_df.reset_index(drop = True)

def col_name_change(col_name,name_dict):
    col_short = re.sub('.*_','',col_name)
    col_idx = re.sub('_.*','',col_name)
    col_idx = int(re.sub('x','',col_idx))
    cat_list = list(name_dict.keys())
    return cat_list[col_idx] + '_' + col_short

def data_clean(df,col_dict_list,description = False):


    pats_df = df[['PATIENT', 'DISEASE', 'ONSET_AGE', 'DEATH_AGE', 'YEARS_TO_DEATH',]].reset_index(drop = True)

    one_hot = df[['MARITAL', 'RACE', 'GENDER']]
    one_hot['MARITAL'] = one_hot['MARITAL'].fillna('S')
    OHE = OneHotEncoder(handle_unknown='ignore')
    hot = OHE.fit_transform(one_hot).toarray()
    hot = pd.DataFrame(hot).rename(columns=lambda x: OHE.get_feature_names()[x])
    name_dict = {i[0]: i[1].unique().tolist() for i in one_hot.iteritems()}
    hot = hot.rename(columns = lambda x: col_name_change(x,name_dict))
    hot_drop_cols = [k + '_' + v[1] for k,v in name_dict.items() if len(v) == 2]
    hot = hot.drop(columns=hot_drop_cols)
    pats_df = pd.concat([pats_df,hot],axis = 1)

    obvs_num_cols = set([j for i in col_dict_list
                     for j in i['obvs_num'] if (find_cols('obvs_num',j,col_dict_list,df) and (len(df[j].dropna())/len(df) >= 0.4) )])

    obvs_num_df = df[obvs_num_cols]
    df_break = [obvs_num_df[df['DISEASE'] == i] for i in df['DISEASE'].unique()]
    imp_list = [imputer(i) for i in df_break]
    obvs_num_df = pd.concat(imp_list).reset_index(drop = True)

    if description:
        med_df = get_cat_df('med_cols',description,col_dict_list,df)
    else:
        med_dict = {}
        for i in range(len(col_dict_list)):
            med_cols = col_dict_list[i]['med_cols']
            sml_med_df = df[med_cols].fillna(0)
            dis_name = df['DISEASE'].unique()[i]
            med_dict[dis_name + '_medication'] = sml_med_df.sum(axis = 1).tolist()
        med_df = pd.DataFrame(med_dict)

    changed_df = [pats_df,med_df,obvs_num_df] + [get_cat_df(i,description,col_dict_list,df) for i in ['obvs_bin', 'cond_cols','proc_cols']]

    full_df = pd.concat(changed_df,axis = 1)

    outcome = (['DEATH_AGE', 'YEARS_TO_DEATH'] +
               list(full_df.filter(regex = 'RATE').columns))



    df_y = full_df[['PATIENT','DISEASE']]
    df_y['code'] = full_df['DISEASE'].copy()
    df_y.loc[:,'code'] = df_y.code.replace(df_y.code.unique(),range(df_y.code.nunique()))
    df_out = full_df[['PATIENT'] + outcome]
    df_X = full_df.drop(['DISEASE'] + outcome, axis = 1)

    return(df_X, df_y, df_out)

def var_count_sorter(var_count_df):
    var_count_df['dif'] = var_count_df.type - var_count_df.fin_row_counts
    if 'MultiIndex' in str(type(var_count_df.index)) :
        dif_indx = var_count_df.loc[var_count_df['dif'] < 0, :].index.tolist()[0]
        var_count_df.index.rename(['noise','var'], inplace = True)
        var_count_df.reset_index(inplace = True)
    else:
        dif_indx = var_count_df.loc[var_count_df['dif'] < 0, ['noise','var']]
        dif_indx = dif_indx.iloc[0].tolist()
    make_up = abs(var_count_df.loc[(var_count_df['noise'] == dif_indx[0]) & (var_count_df['var'] == dif_indx[1])]['dif'].item())
    noise_df = var_count_df[(var_count_df['noise'] == dif_indx[0]) & ~(var_count_df['var'] == dif_indx[1])]
    var_df = var_count_df[~(var_count_df['noise'] == dif_indx[0]) & (var_count_df['var'] == dif_indx[1])]
    else_df = var_count_df[~(var_count_df['noise'] == dif_indx[0]) & ~(var_count_df['var'] == dif_indx[1])]
    df_list = [noise_df,var_df,else_df]
    for i in range(len(df_list)):
        for j in range(len(df_list[i])):
            add = round(make_up/(len(df_list[i]) - j))
            if df_list[i].iloc[j]['dif'] >= add:
                df_list[i].iat[j,5] = df_list[i].iloc[j]['fin_row_counts'] + add
                make_up = make_up - add
            else:
                num = (0 if df_list[i].iloc[j]['dif'] <= 0 else df_list[i].iloc[j]['dif'])
                df_list[i].iat[j,5] = df_list[i].iloc[j]['fin_row_counts'] + num
                make_up = make_up - num
            if make_up <= 0:
                break
    df_list = df_list + [var_count_df[(var_count_df['noise'] == dif_indx[0]) & (var_count_df['var'] == dif_indx[1])]]
    out_df = pd.concat(df_list, sort= False )
    out_df.iat[-1,5] = out_df.iloc[-1]['fin_row_counts'] + out_df.iloc[-1]['dif']
    return(out_df)


def var_ratio_returner(importance,var_n, noise_var_ratio, priority):
    indexes1 = [i for i in ['feature','noise'] for j in noise_var_ratio[0]]
    indexes2 = [j for i in noise_var_ratio for j in ['cont','bin','cat'][:len(noise_var_ratio[0])]]
    var_counter_ex = pd.DataFrame(importance.groupby('noise').type.value_counts())
    var_counter = pd.DataFrame(index = [indexes1,indexes2])
    var_counter = pd.concat([var_counter, var_counter_ex], axis=1).fillna(0)
    var_counter['ratio'] = [i for lists in noise_var_ratio for i in lists]
    var_counter['div'] = round(var_counter.type / var_counter.ratio)

    if var_n == None:
        var_counter['fin_row_counts'] = var_counter.ratio * min(var_counter['div'])
        var_counter.index.rename(['noise', 'var'], inplace=True)
        var_counter.reset_index(inplace=True)
        var_counter_fin = var_counter
    elif priority == 'var_n':
        if var_n > len(importance):
            warnings.warn(
                'There are not enough variables in the data to return the required number of variables (sorry)')
        var_counter['fin_row_counts'] = round(var_counter.ratio * (var_n / sum(var_counter.ratio)))
        count = 0
        while any(var_counter.type < var_counter.fin_row_counts):
            var_counter['dif'] = var_counter.type - var_counter.fin_row_counts
            var_counter = var_count_sorter(var_counter)
            count = count + 1
            print(count)
        warnings.warn('variable ratio could not be maintained while achieving the correct number of vars')
        var_counter_fin = var_counter
    else:
        var_counter['fin_row_counts'] = round(var_counter.ratio * (var_n / sum(var_counter.ratio)))
        var_counter['fin_row_counts'] = list(
            map(
                lambda x: var_counter.iloc[x]['fin_row_counts']if var_counter.iloc[x]['fin_row_counts'] < var_counter.iloc[x]['type'] else var_counter.iloc[x]['type'],
            range(len(var_counter))))
        var_counter.index.rename(['noise', 'var'], inplace=True)
        var_counter.reset_index(inplace=True)
        var_counter_fin = var_counter
    return(var_counter_fin)

def var_getter(var_counter_fin,importance):
    def list_returner(importance,noise,var,n):
        df = importance[(importance['noise'] == noise) & (importance['type'] == var)]
        cols = random.sample(df['vars'].tolist(), k = int(n))
        return(cols)

    var_list = [list_returner(
        importance,
        var_counter_fin.loc[i,'noise'],var_counter_fin.loc[i,'var'],var_counter_fin.loc[i,'fin_row_counts'])
        for i in range(len(var_counter_fin))]
    var_list = [j for i in var_list for j in i]
    return(var_list)


if __name__ == '__main__':

    df_1 = pd.read_csv('test/copd_2021-07-15_4/copdclean.csv')
    df_2 = pd.read_csv('test/dementia_2021-07-15_4/dementiaclean.csv')
    disease_list = [df_1,df_2]
    with open('test/copd_2021-07-15_4/copdcol_dict.json') as f:
        copd_col_dict = json.load(f)

    with open('test/dementia_2021-07-15_4/dementiacol_dict.json') as f:
        dem_col_dict = json.load(f)

    col_dict_list = [copd_col_dict,dem_col_dict]
    df = _combine_disease_dfs(disease_list)
    df_X,df_y,outcomes = data_clean(df,col_dict_list)
