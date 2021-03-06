"""
Functions for generating and cleaning synthea data for one disease
"""

import pandas as pd
import numpy as np
import datetime as dt
import subprocess
import re
import os
import math

def _disease_counter(n,disease,  seed, out_folder = os.getcwd()):
    """
    Generates synthea data for a disease with n number of cases
    :param n: Number of disease cases
    :param disease: synthea disease module as a string
    :param seed: random number for
    :param out_folder: output folder for data (if left blank will go to wd/output/disease)
    :return: creates output in folder (doesn't actually return anything to python) also creates a config file
    showing all the parameters needed to make this so you can save it for later
    """
    # todo write a verbose argument
    # number of patients
    npats = 0

    # conditions not to count

    cvd = ["Coronary Heart Disease",
           "Myocardial Infarction",
           "History of myocardial infarction (situation)",
           "Cardiac Arrest",
           "History of cardiac arrest (situation)",
           "Stroke",
           "Atrial Fibrillation"]
    excepts_dict = {'copd':['Anemia (disorder)'],
                    'dementia':['Pneumonia'],
                    'hypothyrodism':['Anemia']}


    if disease in excepts_dict.keys():
       cvd = cvd + excepts_dict[disease]



    seed_str = str(seed)
    p = "2000"
    conf = "synthea_config.txt"
    count = 0

    # Defining output directory

    file_out = (out_folder+
                '/' +
                re.sub('\*',"",disease, count = 0) +
                "_" +
                str(dt.datetime.now().date()) +
                "_" +
                seed_str
                )


    # run synthea while npats is less then chosen n

    while npats < n and count < 100:
        synth_cmd = ("java -jar synthea-with-dependencies.jar -a 18-100 -m " +
                     disease +
                     " -p " +
                     p +
                     " -s " +
                     seed_str +
                     " -c " +
                     conf +
                     " --exporter.baseDirectory " +
                     file_out)

        run_thing = subprocess.run(synth_cmd , stdout = subprocess.PIPE)
        if run_thing.stdout.decode('utf-8').find('[0 loaded]') > -1:
            raise ValueError('so the disease input doesnt link to a module')

        seed_str = str(int(seed_str) + 1)
        pats = pd.read_csv((file_out + "/csv/conditions.csv"), usecols=['PATIENT','DESCRIPTION'])


        pats = pats[~pats['DESCRIPTION'].isin(cvd)].PATIENT.unique()
        npats_old = npats
        npats = len(pats)



        if npats > 0:
            rate = int(p)/((npats - npats_old) + 1)
        if count > 2 and npats == 0:
            break
        else:
            print(disease, 'count ', count, 'pop ', p, 'pat count ', npats)

        if npats == 0 and int(p) < 200000:
            p = str(int(p) * 5)
        elif npats > 0:
            p = str(min(int(abs(((n - npats) * (rate * 1.05) + 1))), 500000))



        conf = "synthea2_config.txt"
        count = count + 1

    # setting up config file

    set_up = pd.DataFrame({"n": str(n),
              "disease": re.sub('\*', "", disease, count=0),
              "seed": str(seed),
              "count": str(count),
              "pats_n": str(npats)},
                index= [1])

    # outputting config files

    np.savetxt((file_out + "/patstest.txt"),pats, delimiter = ",", fmt = "%s")
    set_up.to_csv((file_out + "/setup.csv"))

def _disease_counter_1d(n,disease,  seed, out_folder = os.getcwd()):
    """
    Generates synthea data for a disease with n number of cases
    :param n: Number of disease cases
    :param disease: synthea disease module as a string
    :param seed: random number for
    :param out_folder: output folder for data (if left blank will go to wd/output/disease)
    :return: creates output in folder (doesn't actually return anything to python) also creates a config file
    showing all the parameters needed to make this so you can save it for later
    """
    # todo write a verbose argument
    # number of patients
    npats = 0

    # conditions not to count

    cvd = ["Coronary Heart Disease",
           "Myocardial Infarction",
           "History of myocardial infarction (situation)",
           "Cardiac Arrest",
           "History of cardiac arrest (situation)",
           "Stroke",
           "Atrial Fibrillation"]
    excepts_dict = {'copd':['Anemia (disorder)'],
                    'dementia':['Pneumonia'],
                    'hypothyrodism':['Anemia']}


    if disease in excepts_dict.keys():
       cvd = cvd + excepts_dict[disease]



    seed_str = str(seed)
    p = "2000"
    conf = "synthea_config.txt"
    count = 0

    # Defining output directory

    file_out = (out_folder+
                '/' +
                re.sub('\*',"",disease, count = 0) +
                "_" +
                str(dt.datetime.now().date()) +
                "_" +
                seed_str
                )


    # run synthea while npats is less then chosen n
    if isinstance(n,int):
        n_dif = [n]
    else:
        n_dif = n
    p_tot = 0
    while any(i > 0 for i in n_dif) and count < 100:
        synth_cmd = ("java -jar synthea-with-dependencies.jar -a 18-100 -m " +
                     disease +
                     " -p " +
                     p +
                     " -s " +
                     seed_str +
                     " -c " +
                     conf +
                     " --exporter.baseDirectory " +
                     file_out)

        run_thing = subprocess.run(synth_cmd , stdout = subprocess.PIPE)
        if run_thing.stdout.decode('utf-8').find('[0 loaded]') > -1:
            raise ValueError('so the disease input doesnt link to a module')

        seed_str = str(int(seed_str) + 1)
        pats = pd.read_csv((file_out + "/csv/conditions.csv"), usecols=['PATIENT','DESCRIPTION'])

        # make sure there are no duplicate rows
        pats = pats[~pats['DESCRIPTION'].isin(cvd)].drop_duplicates()
        #drop pats with more than one disease
        drop_pats = pats.PATIENT.value_counts()
        drop_pats = drop_pats[drop_pats > 1].index.values
        pats = pats[~pats['PATIENT'].isin(drop_pats)]

        rep_n  = pats.DESCRIPTION.nunique()
        if isinstance(n,int):
            n_list = [n for i in range(rep_n)]
        else:
            n_list = n[:rep_n]


        count_list = pd.DataFrame({'npats':pats.DESCRIPTION.value_counts().tolist(), 'n': n_list})
        count_list['ratio'] = count_list['npats']/ count_list['n']
        min_indx = count_list[count_list['ratio'] == count_list['ratio'].min()].index[0]
        npats = count_list.loc[min_indx, 'npats']
        p_tot = p_tot + int(p)
        n_dif = (count_list['n'] - count_list['npats']).to_list()


        if npats > 0:
            rate = p_tot/npats + 1
        if count > 2 and npats == 0:
            break
        else:
            print(disease, 'count ', count, 'pop ', p, 'pat count ',npats)

        if npats == 0 and int(p) < 200000:
            p = str(int(p) * 5)
        elif npats > 0:
            p = str(min(int(abs(min(n_dif) * rate * 1.05 + 1)), 500000))



        conf = "synthea2_config.txt"
        count = count + 1

    # setting up config file

    set_up = pd.DataFrame({"n": str(n),
              "disease": re.sub('\*', "", disease, count=0),
              "seed": str(seed),
              "count": str(count),
              "pats_n": str(npats)},
                index= [1])

    # outputting config files
    def get_labs(df, desc, n):
        df2 = df.copy()
        df2 = df2[df['DESCRIPTION'] == desc][:n]
        return(df2)

    pats_list = [get_labs(pats,pats.DESCRIPTION.unique()[i], n_list[i]) for i in range(len(n_list))]
    pats_list = pd.concat(pats_list, axis = 0, ignore_index= True)
    pats_list.to_csv(file_out + "/patstest.csv")
    set_up.to_csv((file_out + "/setup.csv"))

def _read_files(folder_name, n , description = False, file_list="default", out_file = os.getcwd()):
    """
    Reads files gets all the cases
    :param folder_name: folder where the data lives
    :param file_list: list of files to extract the data, if set to default then it will be the files we think are most
    important (conditions, encounters, medications. observations, patients and procedures)
    :param out_file: folder where to get stuff and put stuff back
    :return: a list of dfs which are all patients, corresponding to file list
    """

    if file_list == "default":
        file_list = ["conditions.csv",
                     "encounters.csv",
                     "medications.csv",
                     "observations.csv",
                     "patients.csv",
                     "procedures.csv"]
    
    # get file locations
    full_file_name = out_file + folder_name + "/csv/"
    get_names = lambda x: re.sub("\.csv", "", x)
    names_list = list(map(get_names, file_list))
    df_list = dict.fromkeys(names_list, 0)
    if description == False:
        patients = pd.read_csv(out_file + folder_name + "/patstest.txt", header=None, names=["PATIENT"]).loc[:(n - 1), :]
    else:
        patients_df = (pd.read_csv(out_file + folder_name + "/patstest.csv").drop(columns = 'Unnamed: 0')
                       .rename({'PATIENT':'PATIENT','DESCRIPTION':'DISEASE'}, axis = 1))
        patients = patients_df.drop(columns = 'DISEASE')
    # sort through files and get only the patients
    for file in file_list:

        # so patients works differently to all the other dfs, this is accounting for that

        if file == "patients.csv":
            df = (pd.read_csv(full_file_name + file, encoding = "latin-1")
                  .pipe(pd.merge, patients, how = "right", left_on = "Id", right_on = "PATIENT"))

        else:
            df = (pd.read_csv(full_file_name + file, encoding = "latin-1")
                  .pipe(pd.merge, patients, how = "inner", on = "PATIENT"))
        file2 = re.sub(".csv", "", file)
        df_list[file2] = df

    if description == False:
        disease = re.sub("_.*", "", folder_name)
        df_list['patients']['DISEASE'] = disease
    else:
         df_list['patients'] = pd.merge(df_list['patients'],patients_df, how = 'outer' ,on = 'PATIENT')

    return (df_list)


# todo export data set
def full_out(disease, df_list, description = False, write_out = False, *args, **kwargs):
    """
    turns the list of data frames into one nice easy to manage data frame
    :param disease: name of disease being looked at
    :param df_list: list of data frames from read files
    :param write_out: whether to write out the data set, set to false or string of folder to write out
    :param args:
    :param kwargs:
    :return: 1 data frame
    """
    def _var_counter(df, onset_df, col_list="Auto", timing=None):
        """
        counts occurance of vars before and after onset
        :param df: data frame
        :param onset_df: data frame containing disease onset
        :param col_list: list of columns to look at
        :param timing: before or after (default none)
        :return: data frame of vars vefore and after disease onset
        """

        if col_list == "Auto":
            col_list = ["PATIENT", "DESCRIPTION"]
        if df.columns.isin(['START']).any():
            df.rename(columns={'START': 'DATE'}, inplace=True)
        df = pd.merge(df, onset_df[['PATIENT', 'START']], how='left', on='PATIENT')

        if timing == 'before':
            df = df[df.DATE <= df.START]
        elif timing == 'after':
            df = df[df.DATE > df.START]
        else:
            timing = ""

        new_df = (df[col_list]
                  .groupby(['PATIENT', 'DESCRIPTION'])
                  .size()
                  .reset_index()
                  .rename(columns={'PATIENT': 'PATIENT', 'DESCRIPTION': 'DESCRIPTION', 0: 'COUNT'})
                  .pivot_table(index='PATIENT', columns='DESCRIPTION', values='COUNT', fill_value=0)
                  .reset_index()
                  .rename(columns=lambda x: x + timing if x != 'PATIENT' else x)
                  )

        return (new_df)

    def _onset_finder(df, except_list="Auto"):
        """
        finds disease onset
        :param df: conditions data frame
        :param except_list: conditions which dont count as diseases
        :return: df with disease onset
        """
        cvd = ["Coronary Heart Disease",
               "Myocardial Infarction",
               "History of myocardial infarction (situation)",
               "Cardiac Arrest",
               "History of cardiac arrest (situation)",
               "Stroke",
               "Atrial Fibrillation"]
        if except_list == "Auto":
            except_list = cvd
        else:
            except_list = except_list + cvd

        df = (df[~df.DESCRIPTION.isin(except_list)]
              .sort_values('START')
              .drop_duplicates('PATIENT'))

        df = df[['PATIENT', 'DESCRIPTION', 'START']]
        return (df)

    def _onset_df_desc(df,pats_df):
        pats = pats_df[['PATIENT','DISEASE']]
        onset_df = pd.merge(df[['PATIENT','DESCRIPTION','START']],
                            pats,
                            how = 'right',
                            left_on = ['PATIENT','DESCRIPTION'],
                            right_on = ['PATIENT','DISEASE']).drop('DISEASE', axis = 1)
        return(onset_df)

    def _year_dif(dates_1, dates_2):
        """
        Takes to pd columns which are dates and finds the difference
        :param dates_1: first date
        :param dates_2: second date (which gets subtracted from the first date)
        :return: column of difference in dates in years
        """
        dif = ((pd.to_datetime(dates_1) - pd.to_datetime(dates_2))
               .astype('str')
               .str.replace(" .*", "")
               .apply(lambda x: float(x) / 365.25 if x != 'NaT' else x))
        return (dif)

    def _obvs_processor(df, onset_df, bin_vars='Auto', vars="All"):
        """
        processes the observations, i knew what it did when I wrote it
        :param df: observations df
        :param onset_df: df with disease onset
        :param bin_vars: variables to binaries
        :param vars: vars to include
        :return: processed observations table with one row per patient
        """
        if vars != "All":
            df = df[df.DESCRIPTION.isin(vars)]

        if bin_vars == 'Auto':
            bin_vars = df.DESCRIPTION[~df.DESCRIPTION.isin(['QALY', 'QOL', 'DALY'])].unique()
        else:
            bin_vars = (bin_vars + df.DESCRIPTION[df.TYPE != 'numeric']).unique()

        df_bin = df[df.DESCRIPTION.isin(bin_vars)]
        df = df[~df.DESCRIPTION.isin(bin_vars)]
        df.loc[:, 'VALUE'] = df.VALUE.astype('float')
        df = (pd.merge(df, onset_df, how='left', on='PATIENT')[['DATE', 'PATIENT', 'DESCRIPTION_x', 'VALUE', 'START']])
        df = df.assign(DATE_DIF=abs(pd.to_datetime(df.DATE) - pd.to_datetime(df.START)))

        df_last = (df.sort_values('DATE', ascending=False)
                   .drop_duplicates(subset=['PATIENT', 'DESCRIPTION_x']))

        df_first = (df.sort_values('DATE_DIF')
                    .drop_duplicates(['PATIENT', 'DESCRIPTION_x']))

        df_last = pd.merge(df_last, df_first, how='outer', on=['PATIENT', 'DESCRIPTION_x'])
        df_last['RATE'] = _year_dif(df_last.DATE_x, df_last.DATE_y)
        df_last['RATE'] = df_last.apply(lambda x: x.VALUE_x - x.VALUE_y / x.RATE if x.RATE != 0 else math.nan, axis=1)
        df_last.apply(lambda x: x.VALUE_x - x.VALUE_y, axis=1)

        df_first = (df_first[['PATIENT', 'DESCRIPTION_x', 'VALUE']]
                    .pivot(index='PATIENT', columns='DESCRIPTION_x', values='VALUE')
                    .rename(columns=lambda x: x + '_FIRST')
                    .reset_index())

        df_last = (df_last[['PATIENT', 'DESCRIPTION_x', 'RATE']]
                   .pivot(index='PATIENT', columns='DESCRIPTION_x', values='RATE')
                   .rename(columns=lambda x: x + '_RATE')
                   .reset_index())

        df_bin = _var_counter(df_bin, onset_df=onset_df)
        df_final = pd.merge(df_first, df_last, how='outer', on='PATIENT')
        df_final = pd.merge(df_final, df_bin, how='outer', on='PATIENT')

        return (df_final)

    def _pats_getter(pats_df, onset_df):
        """
        returns patient df with useful information
        :param pats_df: df with patient information
        :param onset_df: df with onset info
        :return:  patient df with one line per patient and age of onset etc
        """
        df = pd.merge(pats_df, onset_df, how='left', on='PATIENT')
        df = df.assign(ONSET_AGE=_year_dif(df.START, df.BIRTHDATE),
                       DEATH_AGE=_year_dif(df.DEATHDATE, df.BIRTHDATE),
                       YEARS_TO_DEATH=_year_dif(df.DEATHDATE, df.START))
        df = df[['PATIENT', 'DISEASE', 'MARITAL', 'RACE', 'ETHNICITY',
                 'GENDER', 'DESCRIPTION', 'START', 'ONSET_AGE', 'DEATH_AGE',
                 'YEARS_TO_DEATH']]

        return (df)

    if description == False:
        onset_df = _onset_finder(df_list['conditions'])
    else:
        onset_df = _onset_df_desc(df_list['conditions'],df_list['patients'])
    pats_df = _pats_getter(df_list['patients'], onset_df)

    pats_df = pats_df.drop(columns = 'DESCRIPTION')
    obvs_df = _obvs_processor(df_list['observations'], onset_df, bin_vars= '')
    cond_df = df_list['conditions']
    cvd_list = ["Coronary Heart Disease",
                    "Myocardial Infarction",
                    "History of myocardial infarction (situation)",
                    "Cardiac Arrest",
                    "History of cardiac arrest (situation)",
                    "Stroke",
                    "Atrial Fibrillation"]
    excepts_dict = {'copd':['Anemia (disorder)'],
                    'dementia':['Pneumonia'],
                    'hypothyrodism':['Anemia']}


    if disease in excepts_dict.keys():
        disease_list = cvd_list + excepts_dict[disease]
    else:
        disease_list = cvd_list

    df = pd.merge(pats_df, obvs_df, how = 'outer', on = 'PATIENT')
    cond_df = cond_df[cond_df.DESCRIPTION.isin(disease_list)]
    cond_df['DESCRIPTION'] = cond_df.DESCRIPTION.apply(lambda x: 'cvd' if x in cvd_list else x)
    current_cols = df.columns.tolist()
    for i in [cond_df, df_list['procedures'], df_list['encounters'], df_list['medications']]:
        for j in ['before', 'after']:

            x = _var_counter(i, onset_df, timing = j)
            if len(x) > 0:
                df = pd.merge(df, x, how = 'left', on = 'PATIENT')

    new_cols = [i for i in df.columns.tolist() if i not in current_cols]
    df[new_cols] =(df[new_cols]
                   .fillna(0)
                   .apply(lambda x: x.astype('int') , axis = 0))

    if write_out != False:
        df.to_csv(write_out + disease + 'clean.csv')
    return (df)

def _remove_files(path):
    """
    removes the raw synthea data from the folder
    :param path: directory of which to remove stuff
    :return: nothing
    """
    folders = [i for i in os.listdir(path) if '.' not in i]
    [os.remove(path + i + '/' + j) for i in folders for j in os.listdir(path + i)]
    [os.rmdir(path + i) for i in folders]


