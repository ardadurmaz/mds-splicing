#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 17:52:05 2022

@author: ardadurmaz
"""

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import scale
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
import tensorflow as tf


def cust_loss(y_true, y_pred):
    
    y_pred = tf.clip_by_value(y_pred, tf.keras.backend.epsilon(), 1.0-tf.keras.backend.epsilon())
    log_loss = tf.math.add(tf.math.multiply(y_true, tf.math.log(y_pred)), tf.math.multiply(tf.math.multiply(1.0-y_true, tf.math.log(1.0-y_pred)), 0.1))
    
    return tf.math.reduce_sum(tf.math.reduce_mean(-1.0*log_loss, axis=-1), axis=-1)

def cust_loss_mse(y_true, y_pred):
    mse_err = tf.keras.losses.mean_squared_error(y_true, y_pred)
    
    return tf.math.reduce_sum(mse_err, axis=-1)


l_wd = '/media/ardadurmaz/HelperWin/Research/splicing'

# Load Data
psi = scale(np.loadtxt('{}/data/psi_mat.tsv.gz'.format(l_wd), delimiter='\t', dtype=np.float64))
gexp = scale(np.loadtxt('{}/data/expr_mat.tsv.gz'.format(l_wd), delimiter='\t', dtype=np.float64))
mut = np.loadtxt('{}/data/mut_mat.tsv.gz'.format(l_wd), delimiter='\t', dtype=np.float64)

enc = tf.keras.models.load_model('{}/data/splicing_enc_model'.format(l_wd))
dec = tf.keras.models.load_model('{}/data/splicing_dec_model'.format(l_wd))

## Decoder Weights
exp_lat_w = dec.layers[0].get_weights()[0]
exp_lat_b = dec.layers[0].get_weights()[1]

psi_lat_w = dec.layers[2].get_weights()[0]
psi_lat_b = dec.layers[2].get_weights()[1]

mut_lat_w = dec.layers[1].get_weights()[0]
mut_lat_b = dec.layers[1].get_weights()[1]

exp_lat_w_s = dec.layers[3].get_weights()[0]
exp_lat_b_s = dec.layers[3].get_weights()[1]

psi_lat_w_s = dec.layers[5].get_weights()[0]
psi_lat_b_s = dec.layers[5].get_weights()[1]

mut_lat_w_s = dec.layers[4].get_weights()[0]
mut_lat_b_s = dec.layers[4].get_weights()[1]




def local_sig(x=None):
    
    return 1.0/(1.0+np.exp(-1.0*x))



exp_var_res = []
psi_var_res = []
mut_var_res = []


for i in range(3):
    z_mean, z_log_var, z, lat_expr_sec, sd_expr_sec, z_expr, lat_mut_sec, sd_mut_sec, z_mut, lat_psi_sec, sd_psi_sec, z_psi = enc.predict([gexp, mut, psi])
    
    for p in range(z.shape[1]):
        expr_pred = np.dot(z[::,p].reshape([-1,1]), exp_lat_w[p,].reshape([1,-1]))+exp_lat_b
        mut_pred = local_sig(np.dot(z[::,p].reshape([-1,1]), mut_lat_w[p,].reshape([1,-1]))+mut_lat_b)
        mut_pred_null = local_sig(mut_lat_b)
        psi_pred = np.dot(z[::,p].reshape([-1,1]), psi_lat_w[p,].reshape([1,-1]))+psi_lat_b
        
        #mut_pred = local_sig(np.dot(z, mut_lat_w)+mut_lat_b)
            
        expr_ratio = (np.sum(np.var(gexp, axis=0))-np.sum(np.var(gexp-expr_pred, axis=0)))/np.sum(np.var(gexp, axis=0))
        psi_ratio = (np.sum(np.var(psi, axis=0))-np.sum(np.var(psi-psi_pred, axis=0)))/np.sum(np.var(psi, axis=0))
        
        mut_logloss = np.sum((mut*np.log(mut_pred))+((1.0-mut)*np.log((1.0-mut_pred))*0.1))
        mut_logloss_null = np.sum((mut*np.log(mut_pred_null))+((1.0-mut)*np.log((1.0-mut_pred_null))*0.1))
        mut_ratio = 1.0-(mut_logloss/mut_logloss_null)
        
        exp_var_res.append(expr_ratio)
        psi_var_res.append(psi_ratio)
        mut_var_res.append(mut_ratio)
        
        
exp_var_res = np.asarray(exp_var_res)
psi_var_res = np.asarray(psi_var_res)
mut_var_res = np.asarray(mut_var_res)


exp_var_res_all = []
psi_var_res_all = []
mut_var_res_all = []


for i in range(3):
    z_mean, z_log_var, z, lat_expr_sec, sd_expr_sec, z_expr, lat_mut_sec, sd_mut_sec, z_mut, lat_psi_sec, sd_psi_sec, z_psi = enc.predict([gexp, mut, psi])
    
    expr_pred = np.dot(z, exp_lat_w)+exp_lat_b
    mut_pred = local_sig(np.dot(z, mut_lat_w)+mut_lat_b)
    mut_pred_null = local_sig(mut_lat_b)
    psi_pred = np.dot(z, psi_lat_w)+psi_lat_b
    
        
    expr_ratio = (np.sum(np.var(gexp, axis=0))-np.sum(np.var(gexp-expr_pred, axis=0)))/np.sum(np.var(gexp, axis=0))
    psi_ratio = (np.sum(np.var(psi, axis=0))-np.sum(np.var(psi-psi_pred, axis=0)))/np.sum(np.var(psi, axis=0))
    
    mut_logloss = np.sum((mut*np.log(mut_pred))+((1.0-mut)*np.log((1.0-mut_pred))*0.1))
    mut_logloss_null = np.sum((mut*np.log(mut_pred_null))+((1.0-mut)*np.log((1.0-mut_pred_null))*0.1))
    mut_ratio = 1.0-(mut_logloss/mut_logloss_null)
    
    exp_var_res_all.append(expr_ratio)
    psi_var_res_all.append(psi_ratio)
    mut_var_res_all.append(mut_ratio)


exp_var_res_full = []
psi_var_res_full = []
mut_var_res_full = []


for i in range(3):
    z_mean, z_log_var, z, lat_expr_sec, sd_expr_sec, z_expr, lat_mut_sec, sd_mut_sec, z_mut, lat_psi_sec, sd_psi_sec, z_psi = enc.predict([gexp, mut, psi])
    expr_pred, mut_pred, psi_pred = dec.predict([z, z_expr, z_mut, z_psi])    
    mut_pred_null = local_sig(mut_lat_b+mut_lat_b_s)
        
    expr_ratio = (np.sum(np.var(gexp, axis=0))-np.sum(np.var(gexp-expr_pred, axis=0)))/np.sum(np.var(gexp, axis=0))
    psi_ratio = (np.sum(np.var(psi, axis=0))-np.sum(np.var(psi-psi_pred, axis=0)))/np.sum(np.var(psi, axis=0))
    
    mut_logloss = np.sum((mut*np.log(mut_pred))+((1.0-mut)*np.log((1.0-mut_pred))*0.1))
    mut_logloss_null = np.sum((mut*np.log(mut_pred_null))+((1.0-mut)*np.log((1.0-mut_pred_null))*0.1))
    mut_ratio = 1.0-(mut_logloss/mut_logloss_null)
    
    exp_var_res_full.append(expr_ratio)
    psi_var_res_full.append(psi_ratio)
    mut_var_res_full.append(mut_ratio)


exp_var_res_shared = []
psi_var_res_shared = []
mut_var_res_shared = []


for i in range(3):
    z_mean, z_log_var, z, lat_expr_sec, sd_expr_sec, z_expr, lat_mut_sec, sd_mut_sec, z_mut, lat_psi_sec, sd_psi_sec, z_psi = enc.predict([gexp, mut, psi])
    
    expr_pred = np.dot(z_expr, exp_lat_w_s)+exp_lat_b_s
    mut_pred = local_sig(np.dot(z_mut, mut_lat_w_s)+mut_lat_b_s)
    mut_pred_null = local_sig(mut_lat_b_s)
    psi_pred = np.dot(z_psi, psi_lat_w_s)+psi_lat_b_s

        
    expr_ratio = (np.sum(np.var(gexp, axis=0))-np.sum(np.var(gexp-expr_pred, axis=0)))/np.sum(np.var(gexp, axis=0))
    psi_ratio = (np.sum(np.var(psi, axis=0))-np.sum(np.var(psi-psi_pred, axis=0)))/np.sum(np.var(psi, axis=0))
    
    mut_logloss = np.sum((mut*np.log(mut_pred))+((1.0-mut)*np.log((1.0-mut_pred))*0.1))
    mut_logloss_null = np.sum((mut*np.log(mut_pred_null))+((1.0-mut)*np.log((1.0-mut_pred_null))*0.1))
    mut_ratio = 1.0-(mut_logloss/mut_logloss_null)
    
    exp_var_res_full.append(expr_ratio)
    psi_var_res_full.append(psi_ratio)
    mut_var_res_full.append(mut_ratio)


z_mean, z_log_var, z, lat_expr_sec, sd_expr_sec, z_expr, lat_mut_sec, sd_mut_sec, z_mut, lat_psi_sec, sd_psi_sec, z_psi = enc.predict([gexp, mut, psi])
reconstructed = dec.predict([z, z_expr, z_mut, z_psi])

64*np.mean(-1.0*np.mean((mut*np.log(reconstructed[1]))+((1.0-mut)*np.log((1.0-reconstructed[1]))*0.1), axis=1))
64*np.mean(np.mean(np.square(gexp-reconstructed[0]), axis=1))
64*np.mean(np.mean(np.square(psi-reconstructed[2]), axis=1))

mut_w = np.load('{}/data/MutW.npy'.format(l_wd))

exp_var_res_all = []
psi_var_res_all = []
mut_var_res_all = []

for i in range(3):
    z_mean, z_log_var, z, lat_expr_sec, sd_expr_sec, z_expr, lat_mut_sec, sd_mut_sec, z_mut, lat_psi_sec, sd_psi_sec, z_psi = enc.predict([gexp, mut, psi])
    z_expr = np.zeros_like(z_expr)
    z_psi = np.zeros_like(z_psi)
    z_mut = np.zeros_like(z_mut)
    
    reconstructed = dec.predict([z_local, z_expr, z_mut, z_psi])
    mut_pred = np.clip(reconstructed[1], 1.0e-7, 1.0-1.0e-7)
    pred_logloss = -1.0*np.sum((mut*np.log(mut_pred))+((1.0-mut)*np.log((1.0-mut_pred))))
        
    exp_var_res_all.append((np.sum(np.var(gexp, axis=0))-np.sum(np.var(gexp-reconstructed[0], axis=0)))/np.sum(np.var(gexp, axis=0)))
    psi_var_res_all.append((np.sum(np.var(psi, axis=0))-np.sum(np.var(psi-reconstructed[2], axis=0)))/np.sum(np.var(psi, axis=0)))
        
    null_reconstructed = dec.predict([np.zeros_like(z_local), z_expr, z_mut, z_psi])
    null_mut_pred = np.clip(null_reconstructed[1], 1.0e-7, 1.0-1.0e-7)
    null_pred_logloss = -1.0*np.sum((mut*np.log(null_mut_pred))+((1.0-mut)*np.log((1.0-null_mut_pred))))
    mut_var_res_all.append(1.0-(pred_logloss/null_pred_logloss))
        

exp_var_res_all_full = []
psi_var_res_all_full = []
mut_var_res_all_full = []

for i in range(3):
    z_mean, z_log_var, z, lat_expr_sec, sd_expr_sec, z_expr, lat_mut_sec, sd_mut_sec, z_mut, lat_psi_sec, sd_psi_sec, z_psi = enc.predict([gexp, mut, psi])
    
    reconstructed = dec.predict([z_local, z_expr, z_mut, z_psi])
    mut_pred = np.clip(reconstructed[1], 1.0e-7, 1.0-1.0e-7)
    pred_logloss = -1.0*np.sum((mut*np.log(mut_pred)*10.0)+((1.0-mut)*np.log((1.0-mut_pred))))
        
    exp_var_res_all_full.append((np.sum(np.var(gexp, axis=0))-np.sum(np.var(gexp-reconstructed[0], axis=0)))/np.sum(np.var(gexp, axis=0)))
    psi_var_res_all_full.append((np.sum(np.var(psi, axis=0))-np.sum(np.var(psi-reconstructed[2], axis=0)))/np.sum(np.var(psi, axis=0)))
        
    null_reconstructed = dec.predict([np.zeros_like(z_local), z_expr, np.zeros_like(z_mut), z_psi])
    null_mut_pred = np.clip(null_reconstructed[1], 1.0e-7, 1.0-1.0e-7)
    null_pred_logloss = -1.0*np.sum((mut*np.log(null_mut_pred)*10.0)+((1.0-mut)*np.log((1.0-null_mut_pred))))
    mut_var_res_all_full.append(1.0-(pred_logloss/null_pred_logloss))


exp_err_res = np.stack(exp_err_res)
exp_w = dec.layers[0].get_weights()[0]
psi_w = dec.layers[2].get_weights()[0]

#expr_feat = np.unique(np.hstack([np.flip(np.argsort(np.abs(exp_w[i,])))[0:15]for i in range(32)]))
#psi_feat = np.unique(np.hstack([np.flip(np.argsort(np.abs(psi_w[i,])))[0:15]for i in range(32)]))
obs_ft = pd.read_csv('{}/data/combined_annotations.tsv'.format(l_wd), sep='\t')

# Load Data
psi = scale(np.loadtxt('{}/data/psi_mat.tsv.gz'.format(l_wd), delimiter='\t', dtype=np.float32))
gexp = scale(np.loadtxt('{}/data/expr_mat.tsv.gz'.format(l_wd), delimiter='\t', dtype=np.float32))

clust_map = {}
for i in range(15):
    s = 'C-{}'.format(i)
    clust_map[s] = i
    
rf_model = RandomForestClassifier(max_depth=1, class_weight='balanced', n_estimators=1500, oob_score=True, bootstrap=True, random_state=42).fit(psi, [clust_map[s] for s in obs_ft['Cluster']])
local_imp = pd.DataFrame({'importance':rf_model.feature_importances_,
                          'featureIndex':np.arange(psi.shape[1])})
np.savetxt('{}/data/psi_rf_importance.tsv.gz'.format(l_wd), X=local_imp, delimiter='\t')

rf_model = RandomForestClassifier(max_depth=1, class_weight='balanced', n_estimators=1500, oob_score=True, bootstrap=True, random_state=42).fit(gexp, [clust_map[s] for s in obs_ft['Cluster']])
local_imp = pd.DataFrame({'importance':rf_model.feature_importances_,
                          'featureIndex':np.arange(gexp.shape[1])})
np.savetxt('{}/data/expr_rf_importance.tsv.gz'.format(l_wd), X=local_imp, delimiter='\t')