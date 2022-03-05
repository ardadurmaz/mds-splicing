# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 16:15:00 2021

@author: durmaz
"""


from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from sklearn.preprocessing import scale
from sklearn.manifold import MDS
from sklearn.metrics import adjusted_rand_score
from sklearn import mixture
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.mixture import BayesianGaussianMixture
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import OneHotEncoder
from sklearn.manifold import TSNE
import umap
import igraph as ig
import leidenalg as la
import proplot as plot

    

if __name__ == '__main__':    
    
    l_wd = 'Research/splicing/'
    
    # Load Data
    obs_ft = pd.read_csv('{}/data/obs_ft.tsv.gz'.format(l_wd), sep='\t')    
    mut_feat = pd.read_csv('{}/data/mut_mat_vaf_feat.tsv.gz'.format(l_wd), delimiter='\t')
    embedd_res = np.load('{}/data/VAE_Embedding.npy'.format(l_wd))
    
    # Cluster
    bic_res = [] 
    for i in range(30):
        gmm = GaussianMixture(n_components=i+1, covariance_type='diag', n_init=5).fit(embedd_res)
        local_bic = gmm.bic(embedd_res)
        bic_res.append(local_bic)
        
    fig, ax = plot.subplots(nrows=1, ncols=1)
    ax.plot(bic_res)
    ax.format(ylabel='BIC', xlabel='# Clusters')
    fig.savefig('{}/plots/BIC_Plot.pdf'.format(l_wd))
    
    k = np.argmin(bic_res)+1
    clust_labels = GaussianMixture(n_components=k, covariance_type='diag', n_init=10).fit_predict(embedd_res)
    obs_ft['Cluster'] = ['C-{}'.format(i) for i in clust_labels]
    
    obs_ft.to_csv('{}/data/obs_ft_vaf_clusters.tsv.gz'.format(l_wd), sep='\t')
    np.savetxt('{}/data/VAE_Embedding.tsv.gz'.format(l_wd), X=embedd_res, delimiter='\t')

