#!/usr/bin/env python
# coding: utf-8
import pickle
import sys

args = sys.argv

tumor = args[2]
with open(args[1], 'rb') as f:
  obj = pickle.load(f)

print(obj.accuracy_val_list)
