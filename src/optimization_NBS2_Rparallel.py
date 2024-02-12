#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import SRW_vDD as SRW
from multiprocessing.pool import ThreadPool
import itertools
import sys
import pickle
import os

os.chdir("/shares/CIBIO-Storage/BCG/scratch1/ddalfovo/Tools/NBS2_NetPropagation")

args = sys.argv

fold = args[4]
tumor = args[6]
folderSomaticAlt = 'somaticAlterations'

os.mkdir("./results/"+tumor)
os.mkdir("./results/"+tumor+"/pickles/")


edges, features, node_names = SRW.load_network('data/'+folderSomaticAlt+'/'+tumor+'/'+tumor+'_edge2features_PCNet.txt')
P_init_train, sample_names_train = SRW.load_samples('data/'+folderSomaticAlt+'/'+tumor+'/foldOpt/'+tumor+'_alterations_fold'+fold+'_training.tsv', node_names)
P_init_val, sample_names_val = SRW.load_samples('data/'+folderSomaticAlt+'/'+tumor+'/foldOpt/'+tumor+'_alterations_fold'+fold+'_validation.tsv', node_names)
group_labels_train = SRW.load_grouplabels('data/'+folderSomaticAlt+'/'+tumor+'/foldOpt/'+tumor+'_labels_fold'+fold+'_training.tsv')
group_labels_val = SRW.load_grouplabels('data/'+folderSomaticAlt+'/'+tumor+'/foldOpt/'+tumor+'_labels_fold'+fold+'_validation.tsv')

feature_names = []
with open('data/'+folderSomaticAlt+'/'+tumor+'/'+tumor+'_feature_names.txt') as f:
    for line in f.read().rstrip().splitlines():_training
        feature_names.append(line)


feature_names.append('selfloop')
feature_names.append('intercept')

nnodes = len(node_names)

SRW_obj = SRW.SRW_solver(edges, features, nnodes, P_init_train, float(args[1]), group_labels_train, float(args[2]), 
                         w_init_sd=0.01, w=None, feature_names=feature_names, 
                         sample_names=sample_names_train, node_names=node_names, loss='WMW', 
                         norm_type='L1', learning_rate=0.1, update_w_func='Adam', 
                         P_init_val=P_init_val, group_labels_val=group_labels_val, ncpus=10, 
                         maxit=1000, early_stop=500, WMW_b=float(args[3]))
SRW_obj.train_SRW_GD()

with open(args[5]+'/'+tumor+'/pickles/'+'fold'+args[4]+'_'+args[1]+'_'+args[2]+'_'+args[3]+'.pickle', 'wb') as f:
    pickle.dump(SRW_obj, f)
