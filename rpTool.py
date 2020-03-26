#!/usr/bin/env python3

import libsbml
import sys
import logging
from scipy import stats
import numpy as np

sys.path.insert(0, '/home/')
import rpSBML


## Normalise by sigmoidal function
#
#
def nonlin(x,deriv=False):
    if(deriv==True):
        return x*(1-x)
    return 1/(1+np.exp(-x))



## Extract the reaction SMILES from an SBML, query rule_score and write the results back to the SBML
#
# NOTE: all the scores are normalised by their maximal and minimal, and normalised to be higher is better
# Higher is better
#TODO: try to standardize the values instead of normalisation.... Advantage: not bounded
def calculateGlobalScore_json(rpsbml_json,
                              weight_rp_steps,
                              weight_rule_score,
                              weight_fba,
                              weight_thermo,
                              max_rp_steps=15, #TODO: add this as a limit in RP2
                              thermo_ceil=8901.2,
                              thermo_floor=-7570.2,
                              fba_ceil=5.0,
                              fba_floor=0.0,
                              pathway_id='rp_pathway',
                              objective_id='obj_RP1_sink__restricted_biomass',
                              thermo_id='dfG_prime_m'):
    path_norm = {}
    ####################################################################################################### 
    ########################################### REACTIONS #################################################
    ####################################################################################################### 
    for reac_id in rpsbml_json['reactions']:
        for bd_id in rpsbml_json['reactions'][reac_id]:
            ####### Thermo ############
            #lower is better
            #WARNING: we will only take the dfG_prime_m value
            if bd_id[:4]=='dfG_':
                if bd_id not in path_norm:
                    path_norm[bd_id] = []
                try:
                    if thermo_ceil>=rpsbml_json['reactions'][reac_id][bd_id]['value']>=thermo_floor:
                        #min-max feature scaling
                        norm_thermo = (rpsbml_json['reactions'][reac_id][bd_id]['value']-thermo_floor)/(thermo_ceil-thermo_floor)
                    elif rpsbml_json['reactions'][reac_id][bd_id]['value']<thermo_floor:
                        norm_thermo = 0.0
                    elif rpsbml_json['reactions'][reac_id][bd_id]['value']>thermo_ceil:
                        norm_thermo = 1.0
                    norm_thermo = 1.0-norm_thermo
                except (KeyError, TypeError) as e:
                    logging.warning('Cannot find the thermo: '+str(bd_id)+' for the reaction: '+str(reac_id))
                    norm_thermo = 1.0
                rpsbml_json['reactions'][reac_id]['norm_'+bd_id] = {}
                rpsbml_json['reactions'][reac_id]['norm_'+bd_id]['value'] = norm_thermo
                path_norm[bd_id].append(norm_thermo)
            ####### FBA ##############
            #higher is better
            #return all the FBA values
            #------- reactions ----------
            elif bd_id[:4]=='fba_':
                rpsbml_json['fba']['reactions'][reac_id][bd_id] = 0.0
                try:
                    norm_fba = 0.0
                    if fba_ceil>=rpsbml_json['reactions'][reac_id][bd_id]['value']>=fba_floor:
                        #min-max feature scaling
                        norm_fba = (rpsbml_json['reactions'][reac_id][bd_id]['value']-fba_floor)/(fba_ceil-fba_floor)
                    elif rpsbml_json['reactions'][reac_id][bd_id]['value']<=fba_floor:
                        norm_fba = 0.0
                    elif rpsbml_json['reactions'][reac_id][bd_id]['value']>fba_ceil:
                        norm_fba = 1.0
                    rpsbml_json['fba']['reactions'][reac_id][bd_id] = norm_fba
                except (KeyError, TypeError) as e:
                    norm_fba = 0.0
                    logging.warning('Cannot find the objective: '+str(bd_id)+' for the reaction: '+str(reac_id))
                 rpsbml_json['reactions'][reac_id]['norm_'+bd_id] = {}
                 rpsbml_json['reactions'][reac_id]['norm_'+bd_id]['value'] = norm_fba
            elif bd_id=='rule_score':
                if bd_id not in path_norm:
                    path_norm[bd_id] = []
                path_norm[bd_id].append(rpsbml_json['reactions'][reac_id][bd_id]['value'])
            else:
                logging.warning('Not normalising: '+str(bd_id))
    ####################################################################################################### 
    ########################################### PATHWAY ###################################################
    ####################################################################################################### 
    ############### FBA ################
    #higher is better
    for bd_id in rpsbml_json['pathway']:
        if bd_id[:4]=='fba_':
            norm_fba = 0.0
            if fba_ceil>=rpsbml_json['pathway'][bd_id]['value']>=fba_floor:
                #min-max feature scaling
                norm_fba = (rpsbml_json['pathway'][bd_id]['value']-fba_floor)/(fba_ceil-fba_floor)
            elif rpsbml_json['reactions'][reac_id][bd_id]['value']<=fba_floor:
                norm_fba = 0.0
            elif rpsbml_json['reactions'][reac_id][bd_id]['value']>fba_ceil:
                norm_fba = 1.0
            else:
                logging.warning('This flux event should never happen')
            rpsbml_json['pathway']['norm_'+bd_id] = {}
            rpsbml_json['pathway']['norm_'+bd_id]['value'] = norm_fba 
    ############# thermo ################
    #lower is better
    for bd_id in path_norm:
        if bd_id[:4]=='dfG_':
            rpsbml_json['pathway']['norm_'+bd_id] = {}
            rpsbml_json['pathway']['norm_'+bd_id]['value'] = 1.0-sum(path_norm[bd_id])/len(rpsbml_json['reactions'])
    ############# rule score ############
    #higher is better
    for bd_id in path_norm:
        if bd_id=='rule_score':
    if not 'rule_score' in path_norm: 
        logging.warning('Cannot detect rule_score: '+str(path_norm))
        rpsbml_json['pathway']['norm_'+bd_id] = {}
        rpsbml_json['pathway']['norm_'+bd_id]['value'] = 0.0
    else:
        rpsbml_json['pathway']['norm_'+bd_id] = {}
        rpsbml_json['pathway']['norm_'+bd_id]['value'] = sum(path_norm[bd_id])/len(rpsbml_json['reactions'])
    ##### length of pathway ####
    #lower is better
    norm_steps = 0.0
    if len(rpsbml_json['reactions'])>max_rp_steps:
        logging.warning('There are more steps than specified')
        norm_steps = 1.0
    else:
        try:
            norm_steps = (float(len(rpsbml_json['reactions']))-1.0)/(float(max_rp_steps)-1.0)
        except ZeroDivisionError:
            norm_steps = 0.0
    norm_steps = 1.0-norm_steps
    rpsbml_json['pathway']['norm_steps'] = {}
    rpsbml_json['pathway']['norm_steps']['value'] = norm_steps
    ##### global score #########
    try:
        globalScore = (rpsbml_json['pathway']['rule_score']*weight_rule_score+
                       rpsbml_json['pathway']['norm_'+str(thermo_id)]*weight_thermo+
                       rpsbml_json['pathway']['norm_steps']*weight_steps+
                       rpsbml_json['pathway']['norm_'+str(objective_id)]*weight_fba+
                       )/sum([weight_rule_score, weight_thermo, weight_steps, weight_fba])
    except ZeroDivisionError:
        globalScore = 0.0
    except KeyError as e:
        logging.error(e)
        globalScore = 0.0
    rpsbml_json['pathway']['global_score'] = {}
    rpsbml_json['pathway']['global_score']['value'] = globalScore
    return globalScore, rpsbml_json



def genBRSynthJSON(rpsbml, pathway_id='rp_pathway'):
    groups = rpsbml.model.getPlugin('groups')
    fbc = rpsbml.model.getPlugin('fbc')
    rp_pathway = groups.getGroup(pathway_id)
    reactions = rp_pathway.getListOfMembers()
    #Loop through all the reactions
    brsynth_dict = {}
    brsynth_dict['pathway'] = rpsbml.readBRSYNTHAnnotation(rp_pathway.getAnnotation())
    brsynth_dict['reactions'] = {}
    for reaction in reactions:
        reaction = rpsbml.model.getReaction(member.getIdRef())
        annot = reaction.getAnnotation()
        brsynth_dict = rpsbml.readBRSYNTHAnnotation(annot)
        brsynth_dict['reactions'][reaction] = brsynth_dict
    return brsynth_dict


def updateBRSynthPathway(rpsbml, rpsbml_json):
    
