#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 08:08:43 2023

@author: ardadurmaz
"""


import numpy as np
import pandas as pd
from mofapy2.run.entry_point import entry_point
from sklearn.preprocessing import scale


l_wd = '/media/ardadurmaz/HelperWin/Research/splicing'

# Load Data
psi = scale(np.loadtxt('{}/data/psi_mat.tsv.gz'.format(l_wd), delimiter='\t', dtype=np.float32))
gexp = scale(np.loadtxt('{}/data/expr_mat.tsv.gz'.format(l_wd), delimiter='\t', dtype=np.float32))
mut = np.loadtxt('{}/data/mut_mat.tsv.gz'.format(l_wd), delimiter='\t', dtype=np.float32)



ent = entry_point()
data_mat = [[psi], [gexp], [mut]]

ent.set_data_options(scale_groups=False, scale_views=False)
ent.set_data_matrix(data_mat, likelihoods = ["gaussian","gaussian", "bernoulli"])
ent.set_model_options(factors = 50, spikeslab_weights = True, ard_factors = True, ard_weights = True)
ent.set_train_options(iter=1000, convergence_mode="slow", startELBO=1, freqELBO=1, dropR2=0.001, gpu_mode=False, verbose=False, seed=1)

ent.build()
ent.run()
ent.save(outfile='{}/data/mofa_train.hdf5'.format(l_wd))
