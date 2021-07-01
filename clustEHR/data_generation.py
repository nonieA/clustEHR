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


    conditions = {'allergic_rhinitis': ['Perennial allergic rhinitis with seasonal variation',
                           'Perennial allergic rhinitis', 'Seasonal allergic rhinitis'],
     'appendicitis': ['Appendicitis', 'History of appendectomy', 'Rupture of appendix'],
     'asthma': ['Childhood asthma', 'Asthma'],
     'attention_deficit_disorder': ['Child attention deficit disorder'],
     'breast_cancer': ['Malignant neoplasm of breast (disorder)'],
     'bronchitis': ['Acute bronchitis (disorder)'],
     'cerebral_palsy': ['Cerebral palsy (disorder)', 'Spasticity (finding)', 'Epilepsy (disorder)',
                        'Intellectual disability (disorder)', 'Gastroesophageal reflux disease (disorder)',
                        'Poor muscle tone (finding)', 'Excessive salivation (disorder)',
                        'Unable to swallow saliva (finding)', 'Dribbling from mouth (finding)',
                        'Constipation (finding)', 'Pneumonia (disorder)', 'Pain (finding)',
                        'Dislocation of hip joint (disorder)', 'Dystonia (disorder)',
                        'Retention of urine (disorder)'],
     'colorectal_cancer': ['Polyp of colon', 'Primary malignant neoplasm of colon', 'Bleeding from anus',
                           'Protracted diarrhea', 'Secondary malignant neoplasm of colon',
                           'Recurrent rectal polyp', 'Malignant tumor of colon'],
     'congestive_heart_failure': ['Chronic congestive heart failure (disorder)'],
     'copd': ['Pulmonary emphysema (disorder)', 'Chronic obstructive bronchitis (disorder)'],
     'cystic_fibrosis': ['Female Infertility', 'Cystic Fibrosis', 'Diabetes from Cystic Fibrosis',
                         'Infection caused by Staphylococcus aureus',
                         'Sepsis caused by Staphylococcus aureus'],
     'dementia': ["Alzheimer's disease (disorder)","Familial Alzheimer's disease of early onset (disorder)"],
     'dermatitis': ['Contact dermatitis', 'Atopic dermatitis'],
     'ear_infections': ['Otitis media'],
     'epilepsy': ['Seizure disorder', 'History of single seizure (situation)', 'Epilepsy'],
     'fibromyalgia': ['Primary fibromyalgia syndrome'],
     'gallstones': ['Acute Cholecystitis', 'Cholelithiasis'],
     'gout': ['Gout'],
     'hypothyroidism': ['Idiopathic atrophic hypothyroidism'],
     'lung_cancer': ['Suspected lung cancer (situation)', 'Non-small cell lung cancer (disorder)',
                     'Non-small cell carcinoma of lung  TNM stage 1 (disorder)'],
     'opioid_addiction': ['Chronic pain', 'Impacted molars', 'Chronic intractable migraine without aura',
                          'Drug overdose'],
     'osteoarthritis': ['Localized  primary osteoarthritis of the hand', 'Osteoarthritis of hip',
                        'Osteoarthritis of knee'],
     'rheumatoid_arthritis': ['Rheumatoid arthritis'],
     'sinusitis': ['Viral sinusitis (disorder)', 'Chronic sinusitis (disorder)',
                   'Acute bacterial sinusitis (disorder)', 'Sinusitis (disorder)'],
     'sore_throat': ['Acute viral pharyngitis (disorder)', 'Streptococcal sore throat (disorder)'],
     'spina_bifida': ['Spina bifida occulta (disorder)', 'Meningomyelocele (disorder)',
                      'Chiari malformation type II (disorder)', 'Congenital deformity of foot (disorder)'],
     'urinary_tract_infections': ['Recurrent urinary tract infection', 'Cystitis',
                                  'Escherichia coli urinary tract infection']
     }

    condits = conditions[re.sub('\*','',disease)]



    seed_str = str(seed)
    p = "2000"
    conf = "clustEHR\synthea_config.txt"
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
        synth_cmd = ("java -jar clustEHR\synthea-with-dependencies.jar -a 18-100 -p " +
                     p +
                     " -s " +
                     seed_str +
                     " -c " +
                     conf +
                     " --exporter.baseDirectory " +
                     file_out)

        run_thing = subprocess.run(synth_cmd, stdout=subprocess.PIPE)
        if run_thing.stdout.decode('utf-8').find('[0 loaded]') > -1:
            raise ValueError('so the disease input doesnt link to a module')

        seed_str = str(int(seed_str) + 1)
        pats = pd.read_csv((file_out + "/csv/conditions.csv"), usecols=['PATIENT','DESCRIPTION'])


        pats = pats[pats['DESCRIPTION'].isin(condits)].PATIENT.unique()
        npats_old = npats
        npats = len(pats)



        if npats > 0:
            rate = int(p)/((npats - npats_old) + 1)
        if count > 2 and npats == 0:
            break
        else:
            print(disease, 'count ', count, 'pop ', p, 'pat count ', npats)

        if npats == 0 and int(p) < 200000:
            p = str(min(int(p) * 5,10000))
        elif npats > 0:
            p = str(min(int(abs(((n - npats) * (rate * 1.05) + 1))), 10000))



        conf = "clustEHR\synthea2_config.txt"
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
    conditions = {'allergic_rhinitis': ['Perennial allergic rhinitis with seasonal variation',
                           'Perennial allergic rhinitis', 'Seasonal allergic rhinitis'],
     'appendicitis': ['Appendicitis', 'History of appendectomy', 'Rupture of appendix'],
     'asthma': ['Childhood asthma', 'Asthma'],
     'attention_deficit_disorder': ['Child attention deficit disorder'],
     'breast_cancer': ['Malignant neoplasm of breast (disorder)'],
     'bronchitis': ['Acute bronchitis (disorder)'],
     'cerebral_palsy': ['Cerebral palsy (disorder)', 'Spasticity (finding)', 'Epilepsy (disorder)',
                        'Intellectual disability (disorder)', 'Gastroesophageal reflux disease (disorder)',
                        'Poor muscle tone (finding)', 'Excessive salivation (disorder)',
                        'Unable to swallow saliva (finding)', 'Dribbling from mouth (finding)',
                        'Constipation (finding)', 'Pneumonia (disorder)', 'Pain (finding)',
                        'Dislocation of hip joint (disorder)', 'Dystonia (disorder)',
                        'Retention of urine (disorder)'],
     'colorectal_cancer': ['Polyp of colon', 'Primary malignant neoplasm of colon', 'Bleeding from anus',
                           'Protracted diarrhea', 'Secondary malignant neoplasm of colon',
                           'Recurrent rectal polyp', 'Malignant tumor of colon'],
     'congestive_heart_failure': ['Chronic congestive heart failure (disorder)'],
     'copd': ['Pulmonary emphysema (disorder)', 'Chronic obstructive bronchitis (disorder)'],
     'cystic_fibrosis': ['Female Infertility', 'Cystic Fibrosis', 'Diabetes from Cystic Fibrosis',
                         'Infection caused by Staphylococcus aureus',
                         'Sepsis caused by Staphylococcus aureus'],
     'dementia': ["Alzheimer's disease (disorder)","Familial Alzheimer's disease of early onset (disorder)"],
     'dermatitis': ['Contact dermatitis', 'Atopic dermatitis'],
     'ear_infections': ['Otitis media'],
     'epilepsy': ['Seizure disorder', 'History of single seizure (situation)', 'Epilepsy'],
     'fibromyalgia': ['Primary fibromyalgia syndrome'],
     'gallstones': ['Acute Cholecystitis', 'Cholelithiasis'],
     'gout': ['Gout'],
     'hypothyroidism': ['Idiopathic atrophic hypothyroidism'],
     'lung_cancer': ['Suspected lung cancer (situation)', 'Non-small cell lung cancer (disorder)',
                     'Non-small cell carcinoma of lung  TNM stage 1 (disorder)'],
     'opioid_addiction': ['Chronic pain', 'Impacted molars', 'Chronic intractable migraine without aura',
                          'Drug overdose'],
     'osteoarthritis': ['Localized  primary osteoarthritis of the hand', 'Osteoarthritis of hip',
                        'Osteoarthritis of knee'],
     'rheumatoid_arthritis': ['Rheumatoid arthritis'],
     'sinusitis': ['Viral sinusitis (disorder)', 'Chronic sinusitis (disorder)',
                   'Acute bacterial sinusitis (disorder)', 'Sinusitis (disorder)'],
     'sore_throat': ['Acute viral pharyngitis (disorder)', 'Streptococcal sore throat (disorder)'],
     'spina_bifida': ['Spina bifida occulta (disorder)', 'Meningomyelocele (disorder)',
                      'Chiari malformation type II (disorder)', 'Congenital deformity of foot (disorder)'],
     'urinary_tract_infections': ['Recurrent urinary tract infection', 'Cystitis',
                                  'Escherichia coli urinary tract infection']
     }

    condits = conditions[re.sub('\*','',disease)]

    seed_str = str(seed)
    p = "2000"
    conf = "clustEHR\synthea_config.txt"
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
        synth_cmd = ("java -jar clustEHR\synthea-with-dependencies.jar -a 18-100 -p " +
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
        pats = pats[pats['DESCRIPTION'].isin(condits)].drop_duplicates()
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



        conf = "clustEHR\synthea2_config.txt"
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

def cond_finder(cond_dict,var):
    for k,v in cond_dict.items():
        if var in v:
            return(k)

# todo export data set
def full_out(disease, df_list, description = False, write_out = False):
    """
    turns the list of data frames into one nice easy to manage data frame
    :param disease: name of disease being looked at
    :param df_list: list of data frames from read files
    :param write_out: whether to write out the data set, set to false or string of folder to write out
    :param args:
    :param kwargs:
    :return: 1 data frame
    """

    conditions = {
        'allergic_rhinitis': ['Perennial allergic rhinitis with seasonal variation',
                           'Perennial allergic rhinitis', 'Seasonal allergic rhinitis'],
        'appendicitis': ['Appendicitis', 'History of appendectomy', 'Rupture of appendix'],
        'asthma': ['Childhood asthma', 'Asthma'],
        'attention_deficit_disorder': ['Child attention deficit disorder'],
        'breast_cancer': ['Malignant neoplasm of breast (disorder)'],
        'bronchitis': ['Acute bronchitis (disorder)'],
        'cerebral_palsy': ['Cerebral palsy (disorder)', 'Spasticity (finding)', 'Epilepsy (disorder)',
                        'Intellectual disability (disorder)', 'Gastroesophageal reflux disease (disorder)',
                        'Poor muscle tone (finding)', 'Excessive salivation (disorder)',
                        'Unable to swallow saliva (finding)', 'Dribbling from mouth (finding)',
                        'Constipation (finding)', 'Pneumonia (disorder)', 'Pain (finding)',
                        'Dislocation of hip joint (disorder)', 'Dystonia (disorder)',
                        'Retention of urine (disorder)'],
        'colorectal_cancer': ['Polyp of colon', 'Primary malignant neoplasm of colon', 'Bleeding from anus',
                           'Protracted diarrhea', 'Secondary malignant neoplasm of colon',
                           'Recurrent rectal polyp', 'Malignant tumor of colon'],
        'congestive_heart_failure': ['Chronic congestive heart failure (disorder)'],
        'copd': ['Pulmonary emphysema (disorder)', 'Chronic obstructive bronchitis (disorder)'],
        'cystic_fibrosis': ['Female Infertility', 'Cystic Fibrosis', 'Diabetes from Cystic Fibrosis',
                         'Infection caused by Staphylococcus aureus',
                         'Sepsis caused by Staphylococcus aureus'],
        'dementia': ["Alzheimer's disease (disorder)","Familial Alzheimer's disease of early onset (disorder)"],
        'dermatitis': ['Contact dermatitis', 'Atopic dermatitis'],
        'ear_infections': ['Otitis media'],
        'epilepsy': ['Seizure disorder', 'History of single seizure (situation)', 'Epilepsy'],
        'fibromyalgia': ['Primary fibromyalgia syndrome'],
        'gallstones': ['Acute Cholecystitis', 'Cholelithiasis'],
        'gout': ['Gout'],
        'hypothyroidism': ['Idiopathic atrophic hypothyroidism'],
        'lung_cancer': ['Suspected lung cancer (situation)', 'Non-small cell lung cancer (disorder)',
                     'Non-small cell carcinoma of lung  TNM stage 1 (disorder)'],
        'osteoarthritis': ['Localized  primary osteoarthritis of the hand', 'Osteoarthritis of hip',
                        'Osteoarthritis of knee'],
        'rheumatoid_arthritis': ['Rheumatoid arthritis'],
        'sinusitis': ['Viral sinusitis (disorder)', 'Chronic sinusitis (disorder)',
                   'Acute bacterial sinusitis (disorder)', 'Sinusitis (disorder)'],
        'sore_throat': ['Acute viral pharyngitis (disorder)', 'Streptococcal sore throat (disorder)'],
        'spina_bifida': ['Spina bifida occulta (disorder)', 'Meningomyelocele (disorder)',
                      'Chiari malformation type II (disorder)', 'Congenital deformity of foot (disorder)'],
        'urinary_tract_infections': ['Recurrent urinary tract infection', 'Cystitis',
                                  'Escherichia coli urinary tract infection']
    }

    extra_conds = {
        'Obesity': ['Body mass index 30+ - obesity (finding)', 'Body mass index 40+ - severely obese (finding)'],
        'Diabetes': ['Prediabetes', 'Diabetes', 'Metabolic syndrome X (disorder)', 'Hyperglycemia (disorder)',
                     'Diabetic retinopathy associated with type II diabetes mellitus (disorder)',
                     'Diabetic renal disease (disorder)', 'Neuropathy due to type 2 diabetes mellitus (disorder)',
                     'Microalbuminuria due to type 2 diabetes mellitus (disorder)',
                     'Nonproliferative diabetic retinopathy due to type 2 diabetes mellitus (disorder)',
                     'Macular edema and retinopathy due to type 2 diabetes mellitus (disorder)',
                     'Proliferative diabetic retinopathy due to type II diabetes mellitus (disorder)',
                     'Blindness due to type 2 diabetes mellitus (disorder)',
                     'Proteinuria due to type 2 diabetes mellitus (disorder)'],
        'Injury': ['Sprain of ankle', 'Laceration of hand', 'Laceration of forearm', 'First degree burn',
                   'Whiplash injury to neck', 'Facial laceration', 'Bullet wound', 'Laceration of foot',
                   'Laceration of thigh', 'Fracture of clavicle', 'Sprain of wrist', 'Fracture of rib',
                   'Fracture of ankle', 'Fracture subluxation of wrist', 'Tear of meniscus of knee',
                   'Fracture of forearm', 'Pathological fracture due to osteoporosis (disorder)',
                   'Fracture of vertebral column without spinal cord injury', 'Closed fracture of hip',
                   'Second degree burn', 'Injury of tendon of the rotator cuff of shoulder',
                   'Fracture of the vertebral column with spinal cord injury', 'Rupture of patellar tendon',
                   'Third degree burn', 'Injury of anterior cruciate ligament',
                   'Injury of medial collateral ligament of knee', 'Burn injury(morphologic abnormality)',
                   'History of disarticulation at wrist (situation)'],
        'Abnormal pregnancy': ['Miscarriage in first trimester', 'Tubal pregnancy', 'Blighted ovum', 'Preeclampsia',
                               'Fetus with unknown complication', 'Antepartum eclampsia', 'Non-low risk pregnancy'],
        'Head Injury': ['Concussion with no loss of consciousness', 'Concussion with loss of consciousness',
                        'Concussion injury of brain', 'Brain damage - traumatic', 'Traumatic brain injury (disorder)'],
        'Osteoporosis': ['Osteoporosis (disorder)'], 'Hypertension': ['Hypertension'],
        'CVD': ['Hyperlipidemia', 'Atrial Fibrillation', 'Coronary Heart Disease', 'Cardiac Arrest',
                'History of cardiac arrest (situation)', 'Hypertriglyceridemia (disorder)', 'Myocardial Infarction',
                'History of myocardial infarction (situation)'], 'Anemia': ['Anemia (disorder)'],
        'Normal pregnancy': ['Normal pregnancy'],
        'Prostate Cancer': ['Neoplasm of prostate', 'Carcinoma in situ of prostate (disorder)',
                            'Metastasis from malignant tumor of prostate (disorder)'],
        'CKD': ['Chronic kidney disease stage 1 (disorder)', 'Chronic kidney disease stage 2 (disorder)',
                'Chronic kidney disease stage 3 (disorder)'],
        'Substance Abuse': ['Opioid abuse (disorder)', 'Alcoholism','Drug overdose'],
        'Stroke': ['Stroke'],
        'Mental Health Issues': ['Attempted suicide - cut/stab', 'Major depression disorder',
                                 'Posttraumatic stress disorder', 'At risk for suicide (finding)',
                                 'Suicidal deliberate poisoning', 'Major depression  single episode',
                                 'Attempted suicide - suffocation'], 'Allergies': ['Acute allergic reaction'],
        'Paralysis': ['Chronic paralysis due to lesion of spinal cord'], 'Male Infertility': ['Male Infertility'],
        'Sepsis': ['Sepsis caused by Pseudomonas (disorder)']
    }

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

    def _onset_finder(df, condits):
        """
        finds disease onset
        :param df: conditions data frame
        :param except_list: conditions which dont count as diseases
        :return: df with disease onset
        """


        df = (df[df.DESCRIPTION.isin(condits)]
              .sort_values('START')
              .drop_duplicates('PATIENT'))

        df = df[['PATIENT', 'DESCRIPTION', 'START']]
        return(df)

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

    def var_occured(df,var):
        df_edit = df.drop_duplicates(subset = ['PATIENT','DESCRIPTION'])
        desc_list = df_edit['DESCRIPTION'].to_list()
        if desc_list.count(var)/len(df['PATIENT'].unique()) > 0.4:
            return True
        else:
            return False

    def _obvs_processor(df, onset_df, vars="All"):
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


        bin_vars = df.DESCRIPTION[df.TYPE != 'numeric'].unique()
        bin_vars_final = [i for i in bin_vars if var_occured(df,i)]

        non_bin_vars = [i for i in df.DESCRIPTION[df.TYPE == 'numeric'].unique() if var_occured(df,i)]


        df = (pd.merge(df, onset_df, how='left', on='PATIENT')[['DATE', 'PATIENT', 'DESCRIPTION_x', 'VALUE', 'START']])
        df = df.assign(DATE_DIF=abs(pd.to_datetime(df.DATE) - pd.to_datetime(df.START)))
        df_bin = df[df.DESCRIPTION_x.isin(bin_vars_final)]
        df = df[df.DESCRIPTION_x.isin(non_bin_vars)]
        df.loc[:, 'VALUE'] = df.VALUE.astype('float')
        df_last = (df.sort_values('DATE', ascending=False)
                   .drop_duplicates(subset=['PATIENT', 'DESCRIPTION_x']))

        df_first = (df.sort_values('DATE_DIF')
                    .drop_duplicates(['PATIENT', 'DESCRIPTION_x']))

        df_last = pd.merge(df_last, df_first, how='outer', on=['PATIENT', 'DESCRIPTION_x'])
        df_last['RATE'] = _year_dif(df_last.DATE_x, df_last.DATE_y)
        df_last['RATE'] = df_last.apply(lambda x: x.VALUE_x - x.VALUE_y / x.RATE if x.RATE != 0 else math.nan, axis=1)

        df_first = (df_first[['PATIENT', 'DESCRIPTION_x', 'VALUE']]
                    .pivot(index='PATIENT', columns='DESCRIPTION_x', values='VALUE')
                    .rename(columns=lambda x: x + '_FIRST')
                    .reset_index())

        df_last = (df_last[['PATIENT', 'DESCRIPTION_x', 'RATE']]
                   .pivot(index='PATIENT', columns='DESCRIPTION_x', values='RATE')
                   .rename(columns=lambda x: x + '_RATE')
                   .reset_index())

        df_bin = (df_bin.sort_values('DATE_DIF')
                    .drop_duplicates(['PATIENT', 'DESCRIPTION_x']))
        df_bin = df_bin[['PATIENT','VALUE']]
        df_bin['FILL'] = 1
        df_bin = df_bin.pivot_table(index = 'PATIENT',columns='VALUE',values = 'FILL',fill_value = 0).reset_index()

        df_final = pd.merge(df_first, df_last, how='outer', on='PATIENT')
        col_dict = {'obvs_num':list(df_final.drop(columns='PATIENT').columns),
                    'obvs_bin':list(df_bin.drop(columns='PATIENT').columns)}
        df_final = pd.merge(df_final, df_bin, how='outer', on='PATIENT')

        return df_final,col_dict

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
        df = df[['PATIENT', 'DISEASE', 'MARITAL', 'RACE',
                 'GENDER', 'DESCRIPTION', 'START', 'ONSET_AGE', 'DEATH_AGE',
                 'YEARS_TO_DEATH']]

        return (df)

    condits = conditions[disease]

    if description == False:
        onset_df = _onset_finder(df_list['conditions'],condits)
    else:
        onset_df = _onset_df_desc(df_list['conditions'],df_list['patients'])
    pats_df = _pats_getter(df_list['patients'], onset_df)

    pats_df = pats_df.drop(columns = 'DESCRIPTION')
    obvs_df,col_dict = _obvs_processor(df_list['observations'], onset_df)

    df = pd.merge(pats_df, obvs_df, how = 'outer', on = 'PATIENT')

    full_conds_dict = {**extra_conds,**conditions}
    useful_conds = set([j for i in full_conds_dict.values() for j in i])
    cond_df = df_list['conditions']
    cond_df = cond_df[cond_df['DESCRIPTION'].isin(useful_conds)]
    cond_df['DISEASES'] = cond_df['DESCRIPTION'].apply(lambda x:cond_finder(full_conds_dict,x))
    cond_df['VALUE'] = 1
    cond_df = (cond_df[['PATIENT','DISEASES','VALUE']].drop_duplicates(subset = ['PATIENT','DISEASES'])
               .pivot_table(index= 'PATIENT',columns='DISEASES',values='VALUE',fill_value=0)
               .reset_index())

    cond_df = cond_df.drop(columns = disease)
    col_dict['cond_cols'] = list(cond_df.drop(columns='PATIENT').columns)
    proc_df = df_list['procedures']
    proc_df = _var_counter(proc_df,onset_df,timing = '')
    col_dict['proc_cols'] = list(proc_df.drop(columns='PATIENT').columns)
    med_df = df_list['medications']
    med_df = med_df[med_df['REASONDESCRIPTION'].isin(condits +[disease])]
    med_df = _var_counter(med_df,onset_df,timing='')
    col_dict['med_cols'] = list(med_df.drop(columns='PATIENT').columns)
    for i in [cond_df, med_df,proc_df]:
        if len(i) > 0:
            df = pd.merge(df, i, how = 'left', on = 'PATIENT')


    if write_out != False:
        df.to_csv(write_out + disease + 'clean.csv')
    return df,col_dict

def _remove_files(path):
    """
    removes the raw synthea data from the folder
    :param path: directory of which to remove stuff
    :return: nothing
    """
    folders = [i for i in os.listdir(path) if '.' not in i]
    [os.remove(path + i + '/' + j) for i in folders for j in os.listdir(path + i)]
    [os.rmdir(path + i) for i in folders]

if __name__ =='__main__':
    _disease_counter(200,'dementia',4,'trial_data')
    _disease_counter(200,'colorectal_cancer',4,'trial_data')
    _disease_counter_1d(50,'copd',3,'trial_data')
    conf = "clustEHR\\synthea_config.txt"




    file_out = 'trial_data/big_data'


    synth_cmd = ("java -jar clustEHR\synthea-with-dependencies.jar -a 18-100 -p 10000 -s 4 -c " +
                 conf +
                 " --exporter.baseDirectory " +
                 file_out)

    run_thing = subprocess.run(synth_cmd, stdout=subprocess.PIPE)
    if run_thing.stdout.decode('utf-8').find('[0 loaded]') > -1:
        raise ValueError('so the disease input doesnt link to a module')

    df_list = _read_files('dementia_2021-06-28_4',200,out_file='trial_data/')
    df, col_dict = full_out('dementia',df_list)
    df_list1 = _read_files('colorectal_cancer_2021-06-28_4',200,out_file='trial_data/')
    df1,col_dict1 = full_out('colorectal_cancer',df_list1)
    df_list = _read_files('copd_2021-06-28_3',50,description=True,out_file='trial_data/')
    disease = 'copd'
    description = True
    disease_list = [df,df1]
    col_dict_list = [col_dict,col_dict1]
    df, col_dict = full_out('copd',df_list,description=True)
    df_list = [df]
    col_dict_list = [col_dict]