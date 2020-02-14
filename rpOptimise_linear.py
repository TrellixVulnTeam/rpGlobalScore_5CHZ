#!/usr/bin/env python
# coding: utf-8

# In[1]:


import libsbml
import argparse
import os
import io
import tarfile
import glob
import tempfile
import shutil

import sys

import rpSBML
import rpTool
import rpToolServe
import json
import numpy as np
from scipy.optimize import minimize

import rbo


# In[2]:


#print to file and terminal
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("logfile.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass    


# In[3]:


sys.stdout = Logger()


# In[4]:


import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)


# In[5]:


with open('match_scores.json', 'r') as fp:
    match_scores = json.load(fp)


# In[6]:


match_scores


# In[7]:


len_meas_steps = {'measured_2': 9, 'measured_4': 9, 'measured_5': 3, 'measured_6': 4, 'measured_7': 4, 'measured_8': 6, 'measured_9': 4, 'measured_10': 8, 'measured_11': 3, 'measured_13': 4, 'measured_14': 6, 'measured_15': 6, 'measured_17': 4, 'measured_18': 6, 'measured_19': 3, 'measured_20': 3, 'measured_21': 9, 'measured_22': 3, 'measured_23': 2, 'measured_24': 4, 'measured_25': 3, 'measured_26': 6, 'measured_27': 4, 'measured_28': 12, 'measured_29': 2, 'measured_30': 2, 'measured_31': 4, 'measured_32': 3, 'measured_33': 2, 'measured_34': 2, 'measured_35': 4, 'measured_36': 4, 'measured_37': 4, 'measured_38': 4, 'measured_39': 4, 'measured_40': 1, 'measured_41': 1, 'measured_42': 2, 'measured_43': 5, 'measured_44': 6, 'measured_45': 3, 'measured_46': 2, 'measured_47': 4, 'measured_48': 1, 'measured_49': 1, 'measured_50': 1, 'measured_51': 1, 'measured_52': 1, 'measured_53': 2, 'measured_54': 6, 'measured_55': 2, 'measured_56': 2, 'measured_57': 5, 'measured_58': 2, 'measured_59': 2, 'measured_60': 5, 'measured_61': 5, 'measured_62': 3, 'measured_63': 3, 'measured_64': 4, 'measured_65': 4, 'measured_66': 6, 'measured_67': 6, 'measured_68': 3, 'measured_69': 5, 'measured_70': 4, 'measured_71': 2, 'measured_72': 3, 'measured_73': 5, 'measured_74': 4, 'measured_75': 3, 'measured_77': 2, 'measured_78': 2, 'measured_79': 5}

# In[9]:


def matchScore(x):
    all_rbo_ext = []
    for meas in match_scores:
        #retreive the file
        res_tar = glob.glob('/home/mdulac/workspace/Galaxy-SynBioCAD/rpOptimise/results/'+meas+'/rpSelenzyme.tar.xz')[0]
        #compile the list of results
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            top_fileNames = rpToolServe.runGlobalScore_hdd(open(res_tar, 'rb'),
                                                          tmpOutputFolder+'/tmp.tar.xz',
                                                          x[0],
                                                          x[1],
                                                          x[2],
                                                          x[3],
                                                          len_meas_steps[meas],
                                                          10) #this is the topx to return, don't care about it here
            #compile the mean score for the detected pathways
            m_s = {}
            for i in match_scores[meas]:
                m_s[i] = np.mean([match_scores[meas][i]['reac_score'], match_scores[meas][i]['ec_score']])
            rbo_res = rbo.rbo_dict(top_fileNames, m_s, p=0.9)
            all_rbo_ext.append(rbo_res.ext)
            #print('############ '+str(meas)+' ###########')
            #print(top_fileNames)
            #print(m_s)
            #print(rbo_res.ext) #this is the result
    print('############################')
    print(x)
    print(all_rbo_ext)
    print('---> '+str(1.0-np.mean(all_rbo_ext)))
    return 1.0-np.mean(all_rbo_ext)


# In[ ]:


x0 = np.array([1.0, 1.0, 1.0, 1.0])
res = minimize(matchScore, x0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True})

print(res)
