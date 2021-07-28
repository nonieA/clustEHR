import json
import os

from clustEHR.single_data_set_generator import generate_data
import pandas as pd
import numpy as np
import re
import random

def get_noise(noise_dict):
    noise = [list(v.values()) for v in noise_dict.values()]
    return noise

def sort_noise(noise_list,sets):
    first_noise = noise_list[0]
    last_noise = noise_list[-1]
    noise_vals = [get_noise(i) for i in [first_noise,last_noise]]
    if all(k > 0 for i in noise_vals for j in i for k in j):
        first_rat = {i + '_' + 'rat':
                         first_noise['feature'][i]/first_noise['noise'][i] for i in first_noise['feature'].keys()}
        last_rat = {i + '_' + 'rat':
                         last_noise['feature'][i]/last_noise['noise'][i] for i in first_noise['feature'].keys()}

        rat_dict = {
            i:np.linspace(first_rat[i + '_rat'],last_rat[i +'_rat'],num = sets) for i in first_noise['feature'].keys()

        }

        final_dict = [{
            'feature':{
                j:rat_dict[j][i] for j in first_noise['feature'].keys()
            },
            'noise':{
                j:1 for j in first_noise['feature'].keys()
            }
        } for i in range(len(rat_dict['binary']))]

        return final_dict
    else:
        rat_dict = {
            i:{
                j:np.linspace(first_noise[i][j],last_noise[i][j],sets) for j in first_noise[i].keys()
            } for i in first_noise.keys()
        }

        final_dict = [{
            j:{k:rat_dict[j][k][i] for k in first_noise[j].keys()} for j in first_noise.keys()
        } for i in range(sets)]

        return final_dict

def sort_noise_full(noise_list,sets):
    if isinstance(noise_list,list) == False:
        noise_list2 = [noise_list for i in range(sets)]
    elif len(noise_list) == 2:
        noise_list2 = sort_noise(noise_list,sets)
    elif (len(noise_list) != sets) and (noise_list !=2):
        noise_list2 = [noise_list[0] for i in range(sets)]
    else:
        noise_list2 = noise_list
    return [get_noise(i) for i in noise_list2]

def sort_config_list(conf_list,sets):
    if isinstance(conf_list,list) == False:
        return_list = [conf_list for i in range(sets)]
    elif len(conf_list) == sets:
        return_list = conf_list
    elif isinstance(conf_list[0],(int,float)):
        if len(conf_list) == 2:
            return_list = np.linspace(conf_list[0],conf_list[1],sets)
        else:
            return_list = [conf_list[0] for i in range(sets)]
    else:
        return_list = [conf_list[0] for i in range(sets)]
    return return_list

def clean_config(config_file):
    sets = config_file['data_sets']
    new_dict = {k:sort_config_list(v,sets) for k,v in config_file.items() if k not in ['data_sets','noise_var_ratio']}
    new_dict['noise_var_ratio'] = sort_noise_full(config_file['noise_var_ratio'],sets)
    new_dict['clusters'] = [int(i) for i in new_dict['clusters'] if isinstance(i,float)]
    new_dict['seed'] = [int(i) for i in new_dict['seed'] if isinstance(i, float)]
    return new_dict


def clustEHR(config_file):
    out_file = config_file['out_file']
    if os.path.isdir(out_file) == False:
        os.mkdir(out_file)
    config_two = clean_config(config_file)
    config_two['out_file'] = [i + '/data_' + str(idx) for idx,i in enumerate(config_two['out_file'])]
    df = [generate_data(
        n = config_two['n'][i],
        seed = config_two['seed'][i],
        clusters = config_two['clusters'][i],
        vars = config_two['vars'][i],
        noise_var_ratio = config_two['noise_var_ratio'][i],
        var_n= config_two['var_n'][i],
        description= config_two['description'][i],
        priority = config_two['priority'][i],
        out_file = config_two['out_file'][i] ) for i in range(config_file['data_sets'])]
    return df

if __name__ == '__main__':
#    file_list = os.listdir('test')
#    file_list = [i for i in file_list if 'true.json' in i]
#    dict_list = []

#    for i in file_list:
#        with open('test/' +i) as f:
#            new_dict = json.load(f)
#            dict_list.append(new_dict)

#    df_list = [clustEHR(i) for i in dict_list]

    file_list = os.listdir('mc_hammer_test')
    no_noise_dict = {}
    normal_dict = {}

    for i in file_list:
        with open('mc_hammer_test/' +i) as f:
            new_dict = json.load(f)
            name = re.sub('\.json','',i)
            if 'no_noise' in i:
                no_noise_dict[name] = new_dict
            else:
                normal_dict[name] = new_dict

    no_noise_2 = {k + '_2':v for k,v in no_noise_dict.items()}
    for k,v in no_noise_2.items():
        no_noise_2[k]['seed'] = [6,9]

    normal_8 = {k + '_xl':v for k,v in normal_dict.items()}
    for k,v in normal_8.items():
        normal_8[k]['clusters'] = [6,9]
        if 'unequal' not in k:
            normal_8[k]['n'] = [[random.randint(100, 400) for i in range(j)] for j in range(6, 10)]

    normal_dict2 = {**normal_dict,**normal_8}
    normal_dict2 = {k + '_2':v for k,v in normal_dict2.items()}
    for k,v in normal_dict2.items():
        normal_dict2[k]['seed'] = [6,9]

    full_dict = {**no_noise_dict,**no_noise_2,**normal_dict,**normal_8,**normal_dict2}
    for k,v in full_dict.items():
        full_dict[k]['out_file'] = 'mc_hammer_test/' + k

    for k,v in full_dict.items():
        with open('mc_hammer_test/full_configs/' + k +'.json','w') as f:
            json.dump(v,f)