import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import statistics

def rf_noise(df_X,df_y,params = None, og_df = None):
    def rfor_imp(df_X, df_y):
        df_y = df_y.DISEASE
        rfor = RandomForestClassifier(n_estimators=1000).fit(df_X, df_y).feature_importances_
        importance = pd.DataFrame({'vars': df_X.columns, 'importance': rfor})
        importance['rank'] = importance.importance.rank(ascending=False)
        return (importance)

    def list_check(list_1,list_2):
        return(all(elem in list_1 for elem in list_2))


    def de_hot_encode(importance, og_df):
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
    if og_df != None:
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
        if series.unique().tolist() == [0,1]:
            return('bin')
        elif isinstance(series[1],float) or isinstance(series[1],int):
            return('cont')
        else:
            return('cat')


    var_df = pd.DataFrame({'vars':[i for df in args for i in df.columns],
                           'type': [var_bin(df[i]) for df in args for i in df.columns]})
    var_df = var_df.drop_duplicates('vars')
    return(var_df)
