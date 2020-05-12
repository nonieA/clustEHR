import pandas as pd
import os
import random
import data_processing as dp
import re
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LinearRegression
from sklearn import metrics
import statsmodels.api as sm
import seaborn as sns
from plotnine import *
import warnings
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
warnings.simplefilter(action='ignore', category=FutureWarning)

test_data = pd.read_csv('./clustEHR/test3/df_X.csv').drop('Unnamed: 0', axis = 1)
lables = pd.read_csv('./clustEHR/test3/df_y.csv').drop('Unnamed: 0', axis = 1)

# number of clusters
# number of patients
# number of total variables
# number of noise vars

def patient_n(seed,k):
    random.seed(seed)
    n = random.sample([1,k],1)[0]
    pops = random.choices(range(100,501,100),k = n)
    return(pops)

def get_file(file_name):
    disease = re.sub('_[0-9].*','',file_name)
    folder = './clustEHR/test3/' + file_name + '/' + disease + 'clean.csv'
    return(folder)

def vars_get(df_columns,seed):
    random.seed(seed)
    n = random.randint(round(len(df_columns)/2),len(df_columns))
    cols = random.sample(df_columns,n)
    var_n = random.randint(1,round(n/2))
    noise_cols = random.sample(cols,var_n)
    return(cols,noise_cols)

def rfor_imp(df_X, df_y):
    df_y = df_y.code
    df_X = df_X
    rfor = RandomForestClassifier(n_estimators=1000).fit(df_X, df_y).feature_importances_
    importance = pd.DataFrame({'vars': df_X.columns, 'importance': rfor})
    importance['rank'] = importance.importance.rank(ascending=False)
    return (importance)

def rf_import(df,y, noise_cols, seed):
    random.seed(seed)
    df[noise_cols] = df[noise_cols].apply(lambda x: random.sample(x.tolist(),len(x)))
    imp = rfor_imp(df,y)
    imp = imp[imp['vars'].isin(noise_cols)].drop('rank',axis = 1)
    return(imp)

def noise_one(seed):
    random.seed(seed)
    print(seed)
    file_list = os.listdir('./clustEHR/test3')
    file_list = [i for i in file_list if '.csv' not in i]
    clust_n_list = [2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,7,8,8,9,9,10,11,12,13,14,15,16,17,18,19,20]
    clust_n = random.sample(clust_n_list,1)[0]
    disease_files = [get_file(i) for i in random.sample(file_list,clust_n)]
    pats = patient_n(seed,clust_n)
    if len(pats) == 1:
        pats = [pats[0] for i in disease_files]
    disease_list = [pd.read_csv(disease_files[i])[:pats[i]].drop('Unnamed: 0', axis = 1) for i in list(range(clust_n))]
    comb_df, excepts = dp._combine_disease_dfs(disease_list)
    df_X, df_y, outcomes = dp.data_clean(comb_df, comb=excepts)
    all_vars, noise_vars = vars_get(df_X.columns.to_list(),seed)
    df_X = df_X[all_vars]
    imp = [rf_import(df_X,df_y,noise_vars, i) for i in list(range(50))]
    imp2 = pd.concat(imp, axis = 0)
    res = imp2.groupby('vars')['importance'].mean()
    out_dic = {'res': res,
               'dtypes': df_X.dtypes[noise_vars],
               'var_n': len(all_vars),
               'diseases':df_y.DISEASE.unique(),
               'k': clust_n,
               'pats_n':pats}
    return(out_dic)

rf_list = [noise_one(i) for i in list(range(720))]

def change_form(dict_og):
    dict = dict_og.copy()
    if all(i == dict['pats_n'][1] for i in dict['pats_n']):
        dict['pats_n'] = 1
    else:
        dict['pats_n'] = 0
    dis_dicts = {v:1 for v in dict['diseases'].tolist()}
    indx_dict = {'noise_var':dict['res'].index.to_list(),
                 'perc_noise_var':[len(dict['res'])/dict['var_n'] * 100 for i in range(len(dict['res']))]}
    vals_dict = {k:v for k,v in dict.items() if k in ['res', 'dtypes']}
    meta_data_dict = {k:[v for i in range(len(indx_dict['noise_var']))] for k,v in dict.items() if k not in ['res', 'dtypes','diseases']}

    full_dict = {**indx_dict, **vals_dict, **meta_data_dict, **dis_dicts}
    out_df = pd.DataFrame(full_dict)
    return(out_df.reset_index())

rf_list2 = rf_list.copy()

rf_list2 = [change_form(i) for i in rf_list2]

out_df = pd.concat(rf_list2 )
out_df = out_df.drop(['vars','index','noise_var'], axis = 1).fillna(0)

df_y = out_df['res']
df_X = out_df.drop('res', axis = 1)
df_X['dtypes'] = df_X['dtypes'].apply(lambda x: 1 if x == 'int64' else 0)


scaler = MinMaxScaler()
df_X = scaler.fit_transform(df_X)
reg = LinearRegression()
reg = reg.fit(df_X,df_y)
coef = reg.coef_


mod = sm.OLS(df_y,df_X)

fii = mod.fit()

p_values = fii.summary2().tables[1]['P>|t|']


res = pd.DataFrame({'cols': out_df.drop('res', axis = 1).columns.to_list(),'coef': coef,'sig':})

df_X_small = out_df[['dtypes', 'k','pats_n', 'perc_noise_var','var_n']]
df_X_small['dtypes'] = df_X_small['dtypes'].apply(lambda x: 1 if x == 'int64' else 0)
scaler = MinMaxScaler()
df_X_small = scaler.fit_transform(df_X_small)
mod2 = sm.OLS(df_y,df_X_small)
fii = mod2.fit()
coefs = fii.summary2().tables[1]['Coef.']
p_vals = fii.summary2().tables[1]['P>|t|']

res2 = pd.DataFrame({'vars':['dtypes', 'k','pats_n', 'perc_noise_var','var_n'],'coefs':coefs,'p_vals': p_vals})

small_out = out_df[['dtypes', 'k','pats_n', 'var_n']]
small_out.loc[:,'dtypes'] = small_out.loc[:,'dtypes'].apply(lambda x: 1 if x == 'int64' else 0)
scaler = MinMaxScaler()
small_out = scaler.fit_transform(small_out)
df_y = out_df['res']
X_train, X_test, y_train, y_test = train_test_split(small_out, df_y, test_size=0.33, random_state=42)

lm = LinearRegression()
lm = lm.fit(X_train,y_train)
out = lm.predict(X_test)

