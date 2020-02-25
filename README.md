<details><summary>CLICK ME</summary>
<p>
```python
import pandas as pd
clust_data = pd.read_csv('./test_data/output_2020-02-14_16-51-10/cluster_data.csv')
```
</p>
</details>
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

To define the ratio of variable types you use nested lists. 

\[\[number of cont features, number of cat features, number of bin features],
\[number of cont noise variables, number of cat noise variables, number of bin noise variables]]

#### Separation 
You can also control for separation of clusters through using a symmetrical matrix 
of the clusters seperation from each other

\[\[NA, v.close, far],\
\[v.close, NA, medium],\
\[far, medium, NA]]

Disease combinations will then be picked with corresponding seperations. 

## Outputs
There are four CSVs exported: 
#### cluster_data.csv
contains 

## Noise Variables and Variable type 
## Seperation 


