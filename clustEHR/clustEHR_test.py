import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix

def clustering(df,labels):
    k = max(labels) + 1
    x = df.drop(columns = 'PATIENT')
    x_scale = StandardScaler().fit_transform(x)
    x_pca = PCA(n_components=0.9).fit_transform(x_scale)
    km = KMeans(n_clusters=k).fit(x_pca)
    km_labs = km.labels_
    conf = confusion_matrix(labels,km_labs)
    right_list = []
    for i in range(k):
        right_list.append(np.max(conf))
        max_ind = np.argmax(conf)
        row_ind = int((max_ind +1)/k)
        col_ind = max_ind - row_ind*k
        conf[row_ind,:] = 0
        conf[:,col_ind] = 0
    return sum(right_list)/len(labels)

if __name__ == '__main__':
