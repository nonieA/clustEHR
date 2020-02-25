# clustEHR

## About 
This is a package for generating realistic electronic health records specifically for 
testing subtyping methods and evaluation based on data generated using [Synthea EHR 
data generator](https://github.com/synthetichealth/synthea)
## Parameters
#### Number of patients (n): 
This is how many patients per cluster, can be entered as one number (int) or 
a list of numbers corresponding to the number of clusters
#### Clusters: 
To generate clusters with known true labels we have used disease types as a proxy. 
So one cluster will be made up of patients with COPD for example and another will 
have patients with gout. 

These diseases can be specifically selected through inputting a list of disease 
names or randomly selected by puttin in an interger representing hte number of 
clusters 
#### Variables 
You can either select specific variables of define a mix of variable types 
(continous, binary, categorical), the number of variables and the ratio of 
noise variables to features. Sometimes you can not keep the variable type ratio
and have reach the number of desired variables so you can decide which argument 
is more important to preserve.



seperation = None, priority = 'var_n' , 
## Outputs 
## Noise Variables and Variable type 
## Seperation 


