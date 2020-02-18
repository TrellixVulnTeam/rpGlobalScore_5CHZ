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
import sobol_seq
import json

import sys

from scipy import stats

import rpSBML
import rpTool
import rpToolServe
import json
import numpy as np
from scipy.optimize import minimize

import rbo
import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)


all_data_json = {}

'''


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


sys.stdout = Logger()
'''


with open('match_scores.json', 'r') as fp:
    match_scores = json.load(fp)

len_meas_steps = {'measured_2': 9, 'measured_4': 9, 'measured_5': 3, 'measured_6': 4, 'measured_7': 4, 'measured_8': 6, 'measured_9': 4, 'measured_10': 8, 'measured_11': 3, 'measured_13': 4, 'measured_14': 6, 'measured_15': 6, 'measured_17': 4, 'measured_18': 6, 'measured_19': 3, 'measured_20': 3, 'measured_21': 9, 'measured_22': 3, 'measured_23': 2, 'measured_24': 4, 'measured_25': 3, 'measured_26': 6, 'measured_27': 4, 'measured_28': 12, 'measured_29': 2, 'measured_30': 2, 'measured_31': 4, 'measured_32': 3, 'measured_33': 2, 'measured_34': 2, 'measured_35': 4, 'measured_36': 4, 'measured_37': 4, 'measured_38': 4, 'measured_39': 4, 'measured_40': 1, 'measured_41': 1, 'measured_42': 2, 'measured_43': 5, 'measured_44': 6, 'measured_45': 3, 'measured_46': 2, 'measured_47': 4, 'measured_48': 1, 'measured_49': 1, 'measured_50': 1, 'measured_51': 1, 'measured_52': 1, 'measured_53': 2, 'measured_54': 6, 'measured_55': 2, 'measured_56': 2, 'measured_57': 5, 'measured_58': 2, 'measured_59': 2, 'measured_60': 5, 'measured_61': 5, 'measured_62': 3, 'measured_63': 3, 'measured_64': 4, 'measured_65': 4, 'measured_66': 6, 'measured_67': 6, 'measured_68': 3, 'measured_69': 5, 'measured_70': 4, 'measured_71': 2, 'measured_72': 3, 'measured_73': 5, 'measured_74': 4, 'measured_75': 3, 'measured_77': 2, 'measured_78': 2, 'measured_79': 5}

def rpExtractJSONfromSBML(inputTar, pathway_id='rp_pathway'):
    json_out = {}
    with tarfile.open(inputTar, mode='r:xz') as in_tf:
        for member in in_tf.getmembers():
            if not member.name=='':
                member_name = member.name.replace('.xml', '').replace('.rpsbml', '').replace('.sbml', '')
                rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(in_tf.extractfile(member).read().decode('utf-8')))
                json_out[member_name] = {}
                #Groups
                groups = rpsbml.model.getPlugin('groups')
                rp_pathway = groups.getGroup(pathway_id)
                json_out[member_name][pathway_id] = rpsbml.readBRSYNTHAnnotation(rp_pathway.getAnnotation())
                json_out[member_name]['reactions'] = {}
                for reac_mem in rp_pathway.getListOfMembers():
                    reac = rpsbml.model.getReaction(reac_mem.getIdRef())
                    json_out[member_name]['reactions'][reac_mem.getIdRef()] = rpsbml.readBRSYNTHAnnotation(reac.getAnnotation())
                #FBC
                fbc = rpsbml.model.getPlugin('fbc')
                json_out[member_name]['objectives'] = {}
                for objective in fbc.getListOfObjectives():
                   json_out[member_name]['objectives'][objective.getId()] = rpsbml.readBRSYNTHAnnotation(objective.getAnnotation())
    return json_out


def matchScore(x):
    all_rbo = []
    print('############################')
    for meas in match_scores:
        #print('\t'+meas)
        #retreive the file
        res_json = None
        try:
            res_json = all_data_json[meas]
        except KeyError:
            try:
                with open(glob.glob('/home/mdulac/workspace/Galaxy-SynBioCAD/rpOptimise/results/'+meas+'/rpSelenzyme.json')[0]) as json_file:
                    res_json = json.load(json_file)
            except IndexError:
                res_json = rpExtractJSONfromSBML(glob.glob('/home/mdulac/workspace/Galaxy-SynBioCAD/rpOptimise/results/'+meas+'/rpSelenzyme.tar.xz')[0])
                with open('/home/mdulac/workspace/Galaxy-SynBioCAD/rpOptimise/results/'+meas+'/rpSelenzyme.json', 'w') as fp:
                    json.dump(res_json, fp)
            #compile the list of results
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            '''
            top_fileNames = rpToolServe.runGlobalScore_hdd(open(res_tar, 'rb'),
                                                          tmpOutputFolder+'/tmp.tar.xz',
                                                          0.0,
                                                          x[0],
                                                          x[1],
                                                          x[2],
                                                          len_meas_steps[meas],
                                                          10) #this is the topx to return, don't care about it here
            '''
            top_fileNames = rpToolServe.runGlobalScore_json(res_json,
                                                          0.0,
                                                          x[0],
                                                          x[1],
                                                          x[2],
                                                          x[3],
                                                          len_meas_steps[meas],
                                                          10) #this is the topx to return, don't care about it here
            #compile the mean score for the detected pathways
            #print('\t==================== Mean Score ===============')
            #print('\tReaction Match Score: '+str(match_scores[meas]))
            m_s = {}
            for i in match_scores[meas]:
                m_s[i] = np.mean([match_scores[meas][i]['reac_score'], match_scores[meas][i]['ec_score']])
            #print('\tMeasured Match Scores: '+str(m_s))
            #print('\tPredicted Match Scores: '+str(top_fileNames))
            rbo_res = rbo.rbo_dict(top_fileNames, m_s, p=0.9)
            all_rbo.append(rbo_res.ext)
            '''
            all_rbo.append(rbo.RankingSimilarity(
                sorted(top_fileNames, key=top_fileNames.__getitem__), 
                sorted(m_s, key=m_s.__getitem__)).rbo())
            '''
            #print('\tRBO: '+str(rbo_res))
            #print('\t min: '+str(rbo_res.min)+' -- res: '+str(rbo_res.res)+' -- ext: '+str(rbo_res.ext))
            #print('############ '+str(meas)+' ###########')
            #print(top_fileNames)
            #print(m_s)
            #print(rbo_res.ext) #this is the result
    print(x)
    all_opti_res['all_runs'][1.0-np.mean(all_rbo)] = list(x)
    print(all_rbo)
    print('---> '+str(1.0-np.mean(all_rbo)))
    return 1.0-np.mean(all_rbo)


########

print('################################################')
print('################### SOBOL ######################')
print('################################################')

global all_opti_res
all_opti_res = {'all_runs': {}}

#sobol = sobol_seq.i4_sobol_generate(4, 100)
all_sobol = {}
sobol = sobol_seq.i4_sobol_generate(4, 5000)
best_weights = []
best_score = 1.0
for i in sobol:
    #weights = [(y/sum(i))*3.0 for y in i]
    #weights = [y for y in i]
    weights = list(i)
    score = matchScore(weights)
    all_sobol[score] = weights
    if score<best_score:
        best_score = score
        best_weights = weights
    print('-------- BEST ----------')
    print(best_weights)
    print('---> '+str(best_score))

all_opti_res['sobol'] = all_sobol

'''
print('################################################')
print('#################### MINIMISE ##################')
print('################################################')

x0 = np.array(best_weights)
res = minimize(matchScore, x0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True})

#all_opti_res['nelder-mead'] = 
'''

with open('results_sobol.json', 'w') as fp:
    json.dump(all_opti_res, fp)
