# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 09:22:29 2023

@author: Emdadi
"""

from __future__ import division
import os
import numpy as np
import pandas as pd
import scipy.stats
import math

from sklearn import metrics
from sklearn import metrics
from sklearn.metrics import confusion_matrix




def read_inputs(folder):
     
    with open(os.path.join(folder, "Labels_train_test.txt"), "r") as raw:
        
        observation_mat = [line.strip("\n").split() for line in raw]
        
    with open(os.path.join(folder, "similarity_clusters.txt"), "r") as raw:
       # raw.next()
        
        drug_matt = [line.strip("\n").split()[0:] for line in raw]
        
     
    with open(os.path.join(folder, "Pearson_cor_Matrix_train_test.txt"), "r") as raw:
        
        cell_sim_1 = [line.strip("\n").split()[0:] for line in raw]
        
    
       
        
 
    observation_mat = np.array(observation_mat, dtype=np.float64)    
    cell_sim_1 = np.array(cell_sim_1, dtype=np.float64)     
    drug_mat=np.array(drug_matt, dtype=np.float64) 
    drug_mat_main=(drug_mat)
    
   
    
    
    
    return observation_mat, cell_sim_1,cell_sim_1,drug_mat_main


def score_to_exact_rank(s):
   
    p=np.argsort(s)[::-1]
    
    return np.argsort(p)




class LMF:

    def __init__(self, c=1, K1=4, K2=4, r=4, theta=1.3, lambda_p=0.6, lambda_l=0.6, alpha=0.4,beta=0.4, max_iter=50):
        self.c = int(c)  # importance level for positive observations
        self.K1 = int(K1) # number of nearest neighbors used for latent matrix construction
        self.K2 = int(K2) # number of nearest neighbors used for score prediction
        self.r = int(r)  # latent matrices's dimension
        self.theta = float(theta) # Gradien descent learning rate
        self.lambda_p = float(lambda_p) # reciprocal of proteins's variance
        self.lambda_l = float(lambda_l) # reciprocal of subcellulars's variance
        self.alpha = float(alpha) # impact factor of nearest neighbor in constructing method
        self.beta = float(beta)
        self.max_iter = int(max_iter)

    def AGD_optimization(self, seed=None):
        if seed is None:
            self.U = np.sqrt(1/float(self.r))*np.random.normal(size=(self.cells_count, self.r))
            self.V = np.sqrt(1/float(self.r))*np.random.normal(size=(self.locs_count, self.r))
            self.cells_biases = np.sqrt(1/float(self.r))*np.random.normal(size=(self.cells_count, 1))
            self.locs_biases = np.sqrt(1/float(self.r))*np.random.normal(size=(self.locs_count, 1))
        else:
            prng = np.random.RandomState(seed[0])
            prng2= np.random.RandomState(seed[1])
            self.U = prng.normal(loc=0.0, scale=np.sqrt(1/self.lambda_p), size=(self.cells_count, self.r))
            self.V = prng.normal(loc=0.0, scale=np.sqrt(1/self.lambda_l), size=(self.locs_count, self.r))
            self.cells_biases = prng2.normal(size=(self.cells_count, 1))
            self.locs_biases = prng2.normal(size=(self.locs_count, 1))
        cells_sum = np.zeros((self.cells_count, self.U.shape[1]))
        locs_sum = np.zeros((self.locs_count, self.V.shape[1]))
        cells_bias_deriv_sum = np.zeros((self.cells_count, 1))
        locs_bias_deriv_sum = np.zeros((self.locs_count, 1))
        last_log = self.log_likelihood()
        for iter in range(self.max_iter):
            cells,user_bias_deriv = self.deriv(True)
            cells_sum += np.square(cells)
            cells_bias_deriv_sum += np.square(user_bias_deriv)
            vec_step_size = self.theta / (1*np.sqrt(cells_sum))
            bias_step_size = self.theta /(1* np.sqrt(cells_bias_deriv_sum))
            self.cells_biases += bias_step_size * user_bias_deriv
            self.U += vec_step_size * cells
            locs,item_bias_deriv = self.deriv(False)
            locs_sum += np.square(locs)
            locs_bias_deriv_sum += np.square(item_bias_deriv)
            vec_step_size = self.theta / np.sqrt(locs_sum)
            bias_step_size = self.theta /(1* np.sqrt(locs_bias_deriv_sum))
            self.V += vec_step_size * locs
            self.locs_biases += bias_step_size * item_bias_deriv
            curr_log = self.log_likelihood()
            delta_log = (curr_log-last_log)/abs(last_log)
            if abs(delta_log) < 1e-5:
                break
            last_log = curr_log

    def deriv(self, protein):
        if protein:
            vec_deriv = np.dot(self.intMat, self.V)
            bias_deriv = np.expand_dims(np.sum(self.intMat, axis=1), 1)
        else:
            vec_deriv = np.dot(self.intMat.T, self.U)
            bias_deriv = np.expand_dims(np.sum(self.intMat, axis=0), 1)
        SLF = np.dot(self.U, self.V.T) 
        SLF += self.cells_biases
        SLF += self.locs_biases.T
        SLF = np.exp(SLF)
        SLF /= (SLF + self.ones)
        SLF = self.intMat1 * SLF
        if protein:
            vec_deriv -= np.dot(SLF, self.V)
            bias_deriv -= np.expand_dims(np.sum(SLF, axis=1), 1)
            vec_deriv -= self.lambda_p*self.U+self.alpha*np.dot(self.PL, self.U)
        else:
            vec_deriv -= np.dot(SLF.T, self.U)
            bias_deriv -= np.expand_dims(np.sum(SLF, axis=0), 1)
            vec_deriv -= self.lambda_l * self.V+self.beta*np.dot(self.DL, self.V)
        return (vec_deriv, bias_deriv) 

    def log_likelihood(self):
        loglik = 0
        A = np.dot(self.U, self.V.T)
        A += self.cells_biases
        A += self.locs_biases.T
        B = A * self.intMat
        loglik += np.sum(B)
        A=np.array(A, dtype=np.float64)
        A = np.array(np.exp(A), dtype=np.float64)
        A += self.ones

        A = np.log(A)
        A = self.intMat1 * A
        loglik -= np.sum(A)

        loglik -= 0.5 * self.lambda_p * np.sum(np.square(self.U))+0.5 * self.lambda_l * np.sum(np.square(self.V))
        loglik -= 0.5 * self.alpha * np.sum(np.diag((np.dot(self.U.T, self.PL)).dot(self.U)))+0.5 * self.beta * np.sum(np.diag((np.dot(self.V.T, self.DL)).dot(self.V)))
        return loglik

    def construct_neighborhood(self, prots_sim,proteins_sim_2,C_Mat):
        #self.dsMat = C_Mat - np.diag(np.diag(C_Mat))
        self.dsMat = np.array(C_Mat)
        self.PSMat = prots_sim - np.diag(np.diag(prots_sim))
        self.PSMat_2 = proteins_sim_2 - np.diag(np.diag(proteins_sim_2))
        if self.K1 > 0:
            S1 = self.get_nearest_neighbors(self.PSMat, self.K1)
            self.PL = self.laplacian_matrix(S1)
            S2 = self.get_nearest_neighbors(self.dsMat, self.K1)
            self.DL = self.laplacian_matrix(S2)
            S1_2 = self.get_nearest_neighbors(self.PSMat_2, self.K2)
            self.PL_2 = self.laplacian_matrix(S1_2)
            
            
        else:
            self.PL = self.laplacian_matrix(self.PSMat)
            self.DL = self.laplacian_matrix(self.dsMat)
            self.PL_2 = self.laplacian_matrix(self.PSMat_2)
            
            
            
            

    def laplacian_matrix(self, S):
        x = np.sum(S, axis=0)
        y = np.sum(S, axis=1)
       
        L = 0.5*(np.diag(x+y) - (S+S.T))  # neighborhood regularization matrix
        return L

    def get_nearest_neighbors(self, S, size):
        m, n = S.shape
        X = np.zeros((m, n))
        for i in range(m):
            ii = np.argsort(S[i, :])[::-1][:min(size, n)]
            X[i, ii] = S[i, ii]
          #  print(ii)
           # print(len(X[4]))
        return X

    def fix_model(self, W, intMat,C_Mat, cells_sim,cell_lines_sim_2, seed=None):
        
        
        
        self.cells_count, self.locs_count = intMat.shape
        self.ones = np.ones((self.cells_count, self.locs_count))
        self.intMat = self.c*intMat*W
        self.intMat1 = (self.c-1)*intMat * W + self.ones
        x, y = np.where(self.intMat > 0)
        self.train_cells, self.train_locs = set(x.tolist()), set(y.tolist())
        self.construct_neighborhood(cells_sim,cell_lines_sim_2,C_Mat)
        self.AGD_optimization(seed)


    def predict_scores(self, test_data):
        
        prot_ind = np.array(list(self.train_cells))
        train_cells_sim = self.PSMat[:, prot_ind]
        train_cells_sim_2 = self.PSMat_2[:, prot_ind]
        
        scores = []
        for p, l in test_data:
            
            
            ii = np.argsort(train_cells_sim_2[p, :])[::-1][:self.K2]
            
                
            val = np.sum(np.dot(train_cells_sim_2[p, ii], self.U[prot_ind[ii], :])*self.V[l, :])/np.sum(train_cells_sim_2[p, ii])
            val+=np.sum(np.dot(train_cells_sim_2[p, ii], self.cells_biases[prot_ind[ii]]))/np.sum(train_cells_sim_2[p, ii]) + self.locs_biases[l]
            scores.append(np.exp(val)/(1+np.exp(val)))
            
                
        return np.array(scores)
        
    def predict_scores_IC50(self, test_data):
        
        lll_train=[]
        prot_ind = np.array(list(self.train_cells))
        train_cells_sim = self.PSMat[:, prot_ind]
        train_cells_sim_2 = self.PSMat_2[:, prot_ind]
        
        Ic50_values=[]
        p_test=[]
        
        for p, l in test_data:
            ii = np.argsort(train_cells_sim_2[p, :])[::-1][:self.K2]
            val = np.sum(np.dot(train_cells_sim_2[p, ii], self.U[prot_ind[ii], :])*self.V[l, :])/np.sum(train_cells_sim_2[p, ii])
            val+=np.sum(np.dot(train_cells_sim_2[p, ii], self.cells_biases[prot_ind[ii]]))/np.sum(train_cells_sim_2[p, ii]) + self.locs_biases[l]
            Ic50_values.append(val)
           
            
        for i in range (len(test_data)):
            p_test.append(test_data[i])
        
            
          
        
        return  np.array(Ic50_values)    

    def __str__(self):
        return "Model: TCLMF, c:%s, K1:%s, K2:%s, r:%s, lambda_p:%s, lambda_l:%s, alpha:%s, beta:%s,theta:%s, max_iter:%s" % (self.c, self.K1, self.K2, self.r, self.lambda_p, self.lambda_l, self.alpha, self.theta, self.max_iter)
