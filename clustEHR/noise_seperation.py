import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import statistics
import re
import numpy as np
from scipy.stats import mode


def rf_noise(df_X,df_y,params = None, og_df = None):

    def fill_in(ref, col):
        df = pd.DataFrame({'ref': ref.tolist(),
                           'col': col.tolist()})
        maximim = col.max()
        for i in ref.unique():
            if df[df['ref'] == i]['col'].isnull().all():
                df.loc[df['ref'] == i, 'col'] = maximim * -1
            elif df[df['ref'] == i]['col'].isnull().any():
                median = df.loc[df['ref'] == i, 'col'].median()
                df.loc[df['ref'] == i, 'col'].fillna(median)
        return (df['col'].tolist())

    def rfor_imp(df_X, df_y):
        df_y = df_y.code
        df_X = df_X.drop(columns= 'PATIENT')
        rfor = RandomForestClassifier(n_estimators=1000).fit(df_X, df_y).feature_importances_
        importance = pd.DataFrame({'vars': df_X.columns, 'importance': rfor})
        importance['rank'] = importance.importance.rank(ascending=False)
        return (importance)

    def list_check(list_1,list_2):
        return(all(elem in list_1 for elem in list_2))


    def de_hot_encode(importance, og_df = None):
        cat_df = og_df.select_dtypes(include = 'object')
        cat_dict = {i:cat_df[i].unique() for i in cat_df.columns.tolist() if len(cat_df[i].unique()) < len(cat_df[i])}
        x_list = set([re.sub('_.*','',i) for i in importance.vars if re.match('x\d_',i) ])
        var_list = [[[re.sub('.*_','',i) for i in importance.vars if k in i],
                       [importance.loc[importance.vars == i , 'importance'].item() for i in importance.vars if k in i]]
                    for k in x_list]
        cat_dict = {k:[statistics.mean(i[1]) for i in var_list if list_check(v, i[0])] for k,v in cat_dict.items()}
        cat_dict = {k:v[0] for k,v in cat_dict.items() if len(v) > 0}
        cat_df = pd.DataFrame({'vars':list(cat_dict.keys()),'importance': list(cat_dict.values())})
        importance = importance[~importance.vars.str.match('x\d_')]
        importance = pd.concat([importance,cat_df])
        importance['rank'] = importance.importance.rank(ascending=False)
        return(importance)


    importance = rfor_imp(df_X,df_y)
    if og_df is not None:
        importance = de_hot_encode(importance,og_df)

    if params == None:
        params = [0,0.007970]
    if isinstance(params, int):
        params = [0,params]

    importance['noise'] = list(map(lambda x: 'noise' if (x >= params[0] and x <= params[1]) else 'feature',
                                   importance.importance))

    return(importance)

def var_type(*args):
    def var_bin(series):
        if series.unique().tolist() == [0,1] or len(series.unique()) == 2:
            return('bin')
        elif isinstance(series[1],float) or 'int' in str(type(series[1])):
            return('cont')
        else:
            return('cat')

    args = [i for i in args if isinstance(i, pd.core.frame.DataFrame)]
    var_df = pd.DataFrame({'vars':[i for df in args for i in df.columns],
                           'type': [var_bin(df[i]) for df in args for i in df.columns]})
    var_df = var_df.drop_duplicates('vars')
    return(var_df)

def percentile_dif(s1, s2,perc):
    s1_l = np.percentile(s1,perc)
    s1_u = np.percentile(s1,(100 - perc))
    s2_l = np.percentile(s2, perc)
    s2_u = np.percentile(s2, (100 - perc))
    dif = min(abs(s2_u - s1_l),abs(s1_u - s2_l))/max(abs(s2_u - s1_l),abs(s1_u - s2_l))
    return(dif)

def mode_dif(s1, s2):
    if mode(s1) == mode(s2):
        return(0)
    else:
        return(1)




