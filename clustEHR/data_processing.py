"""
Takes multiple disease DFs and processes them to make one clusterable dataset
"""
import pandas as pd
import os
import numpy as np
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.preprocessing import OneHotEncoder
import random
import warnings

def _combine_disease_dfs(disease_list):
    """
    takes individiul diseases and makes them one dataset
    :param disease_list: list of pd.dfs or file paths
    :return: one df and concat dict
    """
    def list_checker(x):
        if isinstance(x, pd.DataFrame):
            return(x)
        elif os.path.isfile(x) and '.csv' in  x:
            return(pd.read_csv(x))
        else:
            raise ValueError('There is a problem in that the input is neither a df or a directory to a df so please can you change that')
    def get_index(col_df):
        drug_indx = [col_df.columns.get_loc(i) for i in col_df.columns if 'before' in i]
        stop = [drug_indx[i] for i in range(1, len(drug_indx)) if drug_indx[i] - drug_indx[(i - 1)] > 1]
        return(stop[-1])

    df_disease_list = [list_checker(i) for i in disease_list]
    int_cols = [df.select_dtypes(include = 'int64').columns.tolist() for df in df_disease_list]
    int_cols = [j for i in int_cols for j in i]

    df = pd.concat(df_disease_list, sort=False)
    df = df.replace([np.inf, -np.inf], np.nan)
    df[int_cols] = (df[int_cols]
                   .fillna(0)
                   .apply(lambda x: x.astype('int'), axis = 0))

    concat_dict = {x.DISEASE[1]:x.iloc[:,get_index(x):].columns.tolist() for x in df_disease_list}
    unique_drug_list = [j for i in concat_dict.values() for j in i]
    unique_drug_list = [i for i in unique_drug_list if unique_drug_list.count(i) == 1]
    def get_unique(drug_list, string,):
        return([i for i in drug_list if i in unique_drug_list and string in i])

    drug_bf = {k + '_drugs_before':get_unique(v,'before') for k,v in concat_dict.items()}
    drug_aft = {k + '_drugs_after':get_unique(v,'after') for k,v in concat_dict.items()}
    concat_dict = {**drug_bf,
                   **drug_aft,
                   'encounters_aft': ['Emergency Encounterafter', 'Encounter for check up (procedure)after',
                                      'Encounter for problemafter', 'Encounter for problem (procedure)after',
                                      'General examination of patient (procedure)after'],
                   'encounters_bf': ['Emergency Encounterafter', 'Encounter for check up (procedure)after',
                                     'Encounter for problemafter', 'Encounter for problem (procedure)after',
                                     'General examination of patient (procedure)after']
                   }

    return(df,concat_dict)


def data_clean(df,
               y_var = ['DISEASE'],
               mult_imps = 'auto',
               comb = None,
               outcome = 'auto',
               zero_rat = 0.9,
               drop_list = ['START', 'ETHNICITY']):


    df_alt = df.copy()
    df_alt = df_alt.drop(drop_list, axis = 1)

    if mult_imps == 'auto':
        mult_imps = ['DALY_FIRST', 'QALY_FIRST', 'QOL_FIRST', 'DALY_RATE', 'QALY_RATE',
                     'QOL_RATE']

    df_imps = df_alt[mult_imps]
    imp = IterativeImputer(max_iter = 100)
    df_imps = imp.fit_transform(df_imps)
    df_imps = pd.DataFrame(df_imps, columns = mult_imps ).reset_index( drop = True)

    df_alt = df_alt.drop(mult_imps, axis = 1)
    df_alt = pd.concat([df_alt.reset_index(drop = True), df_imps],  axis = 1 )

    def column_check(columns, comb_list):
        return([i for i in comb_list if i in columns])

    if comb != None:
        comb = {k:column_check(df_alt.columns,v) for (k,v) in comb.items()}
        comb_dict2 = {k: df[comb.get(k)].sum(axis=1) for k in comb}
        comb_dict2 = pd.DataFrame(comb_dict2).reset_index(drop=True)
        drop_cols = sum(comb.values(), [])

    df_alt = df_alt.drop(columns = drop_cols).join(comb_dict2)



    def zero_sel(df, column, zero_rat, z_type = 'zero'):
        if z_type == 'zero':
            z_count = len(df[df[column] == 0])
        else:
            z_count = df[column].isna().sum()
        if z_count == 0:
            return(False)
        count = len(df)
        if z_count/count > zero_rat:
            return(True)
        else:
            return(False)

    zero_drop = [i for i in df_alt.columns.tolist() if zero_sel(df_alt,i,zero_rat) == True]
    na_drop = [i for i in df_alt.columns.tolist() if zero_sel(df_alt,i,zero_rat, z_type = 'na') == True]
    df_alt = df_alt.drop(zero_drop + na_drop, axis = 1)

    if outcome == 'auto':
        outcome = (['DEATH_AGE', 'YEARS_TO_DEATH'] +
                   list(df_alt.filter(regex = 'aft').columns) +
                   list(df_alt.filter(regex = 'RATE').columns))



    df_y = df_alt[['PATIENT'] + y_var]
    df_y['code'] = df_alt[y_var].copy()
    df_y.loc[:,'code'] = df_y.code.replace(df_y.code.unique(),range(df_y.code.nunique()))
    df_out = df_alt[['PATIENT'] + outcome]
    df_X = df_alt.drop(y_var + outcome, axis = 1)
    df_X['MARITAL'] = df_X['MARITAL'].fillna('S')
    nan_df = [i for i in df_X.columns.tolist() if df_X[i].isna().any() and df_X[i].dtype == 'float' ]
    df_X[nan_df] = df_X[nan_df].apply(lambda x: x.fillna(x.mean()), axis = 0 )
    OHE = OneHotEncoder(handle_unknown='ignore')
    hot = OHE.fit_transform(df_X.select_dtypes(include=['object']).drop(columns = 'PATIENT')).toarray()
    hot = pd.DataFrame(hot).rename(columns=lambda x: OHE.get_feature_names()[x])
    df_X = df_X.reset_index(drop=True)
    df_X = pd.concat([df_X, hot], axis=1)   .select_dtypes(exclude=['object'])

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


