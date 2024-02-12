#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import SRW_vDD_mod as SRW
import pickle

tumor = 'BRCA'
netUsed = 'PCNet'
folderSomaticAlt = 'somaticAlterations'


edges, features, node_names = SRW.load_network('data/'+folderSomaticAlt+'/'+tumor+'/'+tumor+'_edge2features_PCNet.txt')
P_init_train, sample_names_train = SRW.load_samples('data/'+folderSomaticAlt+'/'+tumor+'/'+tumor+'_alterations_training.tsv', node_names)
P_init_val, sample_names_val = SRW.load_samples('data/'+folderSomaticAlt+'/'+tumor+'/'+tumor+'_alterations_validation.tsv', node_names)
group_labels_train = SRW.load_grouplabels('data/'+folderSomaticAlt+'/'+tumor+'/'+tumor+'_labels_training.tsv')
group_labels_val = SRW.load_grouplabels('data/'+folderSomaticAlt+'/'+tumor+'/'+tumor+'_labels_validation.tsv')


feature_names = []
with open('data/'+folderSomaticAlt+'/'+tumor+'/'+tumor+'_feature_names.txt') as f:
    for line in f.read().rstrip().splitlines():
        feature_names.append(line)


feature_names.append('selfloop')
feature_names.append('intercept')
nnodes = len(node_names)

# BRCA
### Put the optimized from 3-fold cross validation step
rst_prob_fix = 0.7 # alpha
lam_fix = 1e-2 # delta
beta_loss_fix = 2e-6 # beta
# PRAD
### Put the optimized from 3-fold cross validation step
rst_prob_fix = 0.1 # alpha
lam_fix = 1 # delta
beta_loss_fix = 2e-3 # beta
# BRCA using 0.5
### Put the optimized from 3-fold cross validation step
rst_prob_fix = 0.3 # alpha
lam_fix = 1e-2 # delta
beta_loss_fix = 2e-5 # beta

SRW_obj = SRW.SRW_solver(edges, features, nnodes, P_init_train, rst_prob_fix, group_labels_train, lam_fix, 
                         w_init_sd=0.01, w=None, feature_names=feature_names, 
                         sample_names=sample_names_train, node_names=node_names, loss='WMW', 
                         norm_type='L1', learning_rate=0.2, update_w_func='Adam', 
                         P_init_val=P_init_val, group_labels_val=group_labels_val, ncpus=50, 
                         maxit=1000, early_stop=500, WMW_b=beta_loss_fix)


SRW_obj.train_SRW_GD()


SRW_obj.generate_Q_and_P_fin()
P_fin_df = SRW_obj.P_fin_df
Q_fin_df = SRW_obj.Q_fin_df

P_fin_df.to_csv('results/'+tumor+'/P_final_mod.txt', sep='\t')
Q_fin_df.to_csv('results/'+tumor+'/Q_final_mod.txt', sep='\t')



with open('results/SRW_obj_'+tumor+'.pickle', 'wb') as f:
    pickle.dump(SRW_obj, f)


# with open('pickle/SRW_obj_.pickle') as f:
    # loaded_obj = pickle.load(f)
