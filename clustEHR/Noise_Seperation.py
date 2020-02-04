def rf_noise(df_X,df_y,noise_vars):
    def rfor_imp(df_X, df_y):
        df_y = df_y.DISEASE
        rfor = RandomForestClassifier(n_estimators=1000).fit(df_X, df_y).feature_importances_
        importance = pd.DataFrame({'var': df_X.columns, 'importance': rfor})
        importance['rank'] = importance.importance.rank(ascending=False)
        return (importance)
    rf_imp = rfor_imp(df_X,df_y)

    if params == None:
        params = [0,0.007970]
    if isinstance(params, int):
        params = [0,params]

    importance['noise'] = list(map(lambda x: 'noise' if x >= params[0] and <= params[1] else 'feature' ))

    return(importance)

def var_type(df):
    def var_bin(series):
        if series.unique().tolist() == [0,1]:
            return('bin')
        elif isinstance(series,float) or isinstance(series,int):
            return('cont')
        else:
            return('cat')
    var_df = pd.DataFrame({'var':df.columns.tolist(), 'type': [var_bin(df[i]) for i in df.columns]})
    return(var_df)

def