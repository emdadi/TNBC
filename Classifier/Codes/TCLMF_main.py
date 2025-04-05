# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 09:20:46 2023

@author: Emdadi
"""

from __future__ import division
import os
import numpy as np
import pandas as pd


from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import  RepeatedKFold
from TCLMF_functions import read_inputs,score_to_exact_rank
from TCLMF_functions import LMF

data_folder = os.path.join(os.path.pardir, 'Datasets') 

observation_mat, proteins_sim,proteins_sim_2,drugMat=read_inputs(data_folder) 
 
seed = [80162,45929]


model = LMF(c=10, K1=7, K2=7, r=40, lambda_p=0.7, lambda_l=0.2, alpha=1.0,beta=0.0, theta=1.5, max_iter=1000)

   
test_indexx=[]
train_indexx=[]
test_indexx.append(116) 
test_index=np.array(test_indexx) 


for i in range(116):
     train_indexx.append(i)   
 
train_index=np.array(train_indexx)  


   
test_location_mat = np.array(observation_mat)
   
test_location_mat[train_index] = 0
train_location_mat = np.array(observation_mat - test_location_mat)
   
true_result = np.array(test_location_mat[test_index])
   
    
x = np.repeat(test_index, len(observation_mat[0]))
y = np.arange(len(observation_mat[0]))
  
y = np.tile(y, len(test_index))
  
model.fix_model(train_location_mat, train_location_mat,drugMat, proteins_sim,proteins_sim_2, seed)
    
scores = np.reshape(model.predict_scores(zip(x, y)), true_result.shape)
#pred_Ic50_values=np.reshape(model.predict_scores_IC50(zip(x, y)), true_result.shape)
   
print('probability',(scores[0]))
#print('Ic50_values',(pred_Ic50_values[0]))