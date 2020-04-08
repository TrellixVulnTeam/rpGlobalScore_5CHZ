#!/usr/bin/env python3

import libsbml
import sys
import logging
from scipy import stats
import numpy as np

import rpSBML


## Normalise by sigmoidal function
#
#
def nonlin(x,deriv=False):
    if(deriv==True):
        return x*(1-x)
    return 1/(1+np.exp(-x))

#TODO: use rpSBML annottions function instead of manual
## Extract the reaction SMILES from an SBML, query selenzyme and write the results back to the SBML
#
# NOTE: all the scores are normalised by their maximal and minimal, and normalised to be higher is better
# Higher is better
#TODO: try to standardize the values instead of normalisation.... Advantage: not bounded
def calculateGlobalScore(rpsbml,
                         weight_rp_steps,
                         weight_rule_score,
                         weight_fba,
                         weight_thermo,
                         max_rp_steps, #fix this to 15 or something
                         thermo_ceil=8901.2,
                         thermo_floor=-7570.2,
                         fba_ceil=5.0,
                         fba_floor=0.0,
                         pathway_id='rp_pathway',
                         objective_id='obj_RP1_sink__restricted_biomass',
                         thermo_id='dfG_prime_m'):
    groups = rpsbml.model.getPlugin('groups')
    fbc = rpsbml.model.getPlugin('fbc')
    rp_pathway = groups.getGroup(pathway_id)
    members = rp_pathway.getListOfMembers()
    reactions_data = {'rule_score': {'reactions': {}, 'global': 0.0},
                      'fba': {'reactions': {}, 'global': {}},
                      'thermo': {'reactions': {}, 'global': {}}}
    #Loop through all the reactions
    for member in members:
        reactions_data['rule_score']['reactions'][member.getIdRef()] = 0.0
        reactions_data['fba']['reactions'][member.getIdRef()] = {}
        reactions_data['thermo']['reactions'][member.getIdRef()] = {}
        reaction = rpsbml.model.getReaction(member.getIdRef())
        annot = reaction.getAnnotation()
        brsynth_dict = rpsbml.readBRSYNTHAnnotation(annot)
        ####### rule_score ###########
        #higher is better
        reactions_data['rule_score']['reactions'][member.getIdRef()] = float(brsynth_dict['rule_score']['value'])
        reactions_data['rule_score']['global'] += float(brsynth_dict['rule_score']['value'])
        ####### Thermo ############
        #lower is better
        #WARNING: we will only take the dfG_prime_m value
        for bd_id in brsynth_dict:
            if bd_id[:4]=='dfG_':
                reactions_data['thermo']['reactions'][member.getIdRef()][bd_id] = 0.0
                if not bd_id in reactions_data['thermo']['global'].keys():
                    reactions_data['thermo']['global'][bd_id] = 0.0
                #reactions_data['thermo']['global'][bd_id] = 0.0
                try:
                    if thermo_ceil>=brsynth_dict[bd_id]['value']>=thermo_floor:
                        norm_thermo = (brsynth_dict[bd_id]['value']-thermo_floor)/(thermo_ceil-thermo_floor)
                    elif brsynth_dict[bd_id]['value']<thermo_floor:
                        norm_thermo = 0.0
                    elif brsynth_dict[bd_id]['value']>thermo_ceil:
                        norm_thermo = 1.0
                    norm_thermo = 1.0-norm_thermo
                    reactions_data['thermo']['reactions'][member.getIdRef()][bd_id] = norm_thermo
                    reactions_data['thermo']['global'][bd_id] += norm_thermo
                    rpsbml.addUpdateBRSynth(reaction, 'norm_'+bd_id, norm_thermo)
                except (KeyError, TypeError) as e:
                    logging.warning('Cannot find the thermo: '+str(bd_id)+' for the reaction: '+str(member.getIdRef()))
                    rpsbml.addUpdateBRSynth(reaction, 'norm_'+bd_id, 0.0)
        ####### FBA ##############
        #higher is better
        #return all the FBA values
        norm_fba = 0.0
        for bd_id in brsynth_dict:
            if bd_id[:4]=='fba_':
                reactions_data['fba']['reactions'][member.getIdRef()][bd_id] = 0.0
                try:
                    if fba_ceil>=brsynth_dict[bd_id]['value']>=fba_floor:
                        #min-max feature scaling
                        norm_fba = (brsynth_dict[bd_id]['value']-fba_floor)/(fba_ceil-fba_floor)
                    elif brsynth_dict[bd_id]['value']<=fba_floor:
                        norm_fba = 0.0
                    elif brsynth_dict[bd_id]['value']>fba_ceil:
                        norm_fba = 1.0
                    reactions_data['fba']['reactions'][member.getIdRef()][bd_id] = norm_fba
                    rpsbml.addUpdateBRSynth(reaction, 'norm_'+bd_id, norm_fba)
                except (KeyError, TypeError) as e:
                    logging.warning('Cannot find the objective: '+str(bd_id)+' for the reaction: '+str(member.getIdRef()))
                    rpsbml.addUpdateBRSynth(reaction, 'norm_'+bd_id, 0.0)
    ##############################
    ##### target FBA value #######
    ##############################
    #higher is better
    #loop through all the different objectives and normalise the values
    #find the objective
    for objective in fbc.getListOfObjectives():
        brsynth_dict = rpsbml.readBRSYNTHAnnotation(objective.getAnnotation())
        if not objective.getId() in reactions_data['fba']['global'].keys():
            reactions_data['fba']['global'][objective.getId()] = 0.0
        try:
            if fba_ceil>=float(brsynth_dict['flux_value']['value'])>=fba_floor:
                norm_fba = (round(float(brsynth_dict['flux_value']['value']), 4)-fba_floor)/(fba_ceil-fba_floor)
            elif float(brsynth_dict['flux_value']['value'])<fba_floor:
                norm_fba = 0.0
            elif float(brsynth_dict['flux_value']['value'])>fba_ceil:
                norm_fba = 1.0
            reactions_data['fba']['global'][objective.getId()] = norm_fba
        except (KeyError, TypeError) as e:
            logging.warning('Could not retreive flux value: '+str(objective.getId()))
            continue
        rpsbml.addUpdateBRSynth(objective, 'norm_flux_value', norm_fba)
        rpsbml.addUpdateBRSynth(rp_pathway, 'norm_'+objective.getId(), norm_fba)
        for flux_obj in objective.getListOfFluxObjectives():
            brsynth_dict = rpsbml.readBRSYNTHAnnotation(flux_obj.getAnnotation())
            try:
                norm_fba = (round(float(brsynth_dict['flux_value']['value']), 4)-fba_floor)/(fba_ceil-fba_floor)
            except (KeyError, TypeError) as e:
                norm_fba = 0.0
            rpsbml.addUpdateBRSynth(flux_obj, 'norm_flux_value', norm_fba)
    try:
        target_norm_fba = reactions_data['fba']['global'][objective_id]
    except KeyError:
        logging.warning('Detected a key error for '+str(objective_id)+' in target_norm_fba calculation')
        target_norm_fba = 0.0
    brsynth_dict = rpsbml.readBRSYNTHAnnotation(rp_pathway.getAnnotation())
    for bd_id in brsynth_dict:
        if bd_id[:4]=='dfG_':
            try:
                reactions_data['thermo']['global'][bd_id] = reactions_data['thermo']['global'][bd_id]/float(len(members))
                rpsbml.addUpdateBRSynth(rp_pathway, 'norm_'+bd_id, reactions_data['thermo']['global'][bd_id])
            except (KeyError, TypeError) as e:
                rpsbml.addUpdateBRSynth(rp_pathway, 'norm_'+bd_id, 0.0)
                logging.warning('Cannot find the thermo: '+str(bd_id)+' for the reaction: '+str(bd_id))
    try:
        target_norm_thermo_mean = np.mean([reactions_data['thermo']['reactions'][i][thermo_id] for i in reactions_data['thermo']['reactions']])
    except KeyError:
        target_norm_thermo_mean = 0.0
    ############################
    ##### length of members ####
    ############################
    #lower is better
    if len(members)>max_rp_steps:
        logging.warning('There are more steps than specified')
        norm_steps = 1.0
    else:
        try:
            norm_steps = float(len(members))/float(max_rp_steps)
        except ZeroDivisionError:
            norm_steps = 0.0
    norm_steps = 1.0-norm_steps
    rpsbml.addUpdateBRSynth(rp_pathway, 'norm_steps', norm_steps)
    ############################
    ####### rule_score #########
    ############################
    #higher is better
    #WARNING: using norm here is redundant but is used to detect the type of annotation
    norm_rScore = float(len(members))/reactions_data['rule_score']['global']
    rpsbml.addUpdateBRSynth(rp_pathway, 'norm_rule_scores', norm_rScore)
    ############################
    ##### global score #########
    ############################
    try:
        globalScore = (norm_steps*weight_rp_steps+
                       target_norm_fba*weight_fba+
                       norm_rScore*weight_rule_score+
                       target_norm_thermo_mean*weight_thermo
                       )/sum([weight_rp_steps, weight_fba, weight_thermo, weight_rule_score])
    except ZeroDivisionError:
        logging.warning('Global score calculation is dividing everything by 0')
        globalScore = 0.0
    rpsbml.addUpdateBRSynth(rp_pathway, 'global_score', globalScore)
    return globalScore



#TODO: use rpSBML annottions function instead of manual
## Extract the reaction SMILES from an SBML, query selenzyme and write the results back to the SBML
#
# NOTE: all the scores are normalised by their maximal and minimal, and normalised to be higher is better
# Higher is better
#TODO: try to standardize the values instead of normalisation.... Advantage: not bounded
def calculateGlobalScore_json(rpsbml_json,
                              weight_rp_steps,
                              weight_rule_score,
                              weight_fba,
                              weight_thermo,
                              max_rp_steps, #fix this to 15 or something
                              thermo_ceil=8901.2,
                              thermo_floor=-7570.2,
                              fba_ceil=5.0,
                              fba_floor=0.0,
                              pathway_id='rp_pathway',
                              objective_id='obj_RP1_sink__restricted_biomass',
                              thermo_id='dfG_prime_m'):
    reactions_data = {'rule_score': {'reactions': {}, 'global': 0.0},
                      'fba': {'reactions': {}, 'global': {}},
                      'thermo': {'reactions': {}, 'global': {}}}
    members = len(rpsbml_json['reactions'])
    #Loop through all the reactions
    #reactions
    for reac_id in rpsbml_json['reactions']:
        reactions_data['rule_score']['reactions'][reac_id] = 0.0
        reactions_data['fba']['reactions'][reac_id] = {}
        reactions_data['thermo']['reactions'][reac_id] = {}
        brsynth_dict = rpsbml_json['reactions'][reac_id]
        ####### rule_score ###########
        #higher is better
        reactions_data['rule_score']['reactions'][member.getIdRef()] = float(brsynth_dict['rule_score'])
        reactions_data['rule_score']['global'] += float(brsynth_dict['rule_score'])
        ####### Thermo ############
        #lower is better
        #WARNING: we will only take the dfG_prime_m value
        for bd_id in brsynth_dict:
            if bd_id[:4]=='dfG_':
                reactions_data['thermo']['reactions'][reac_id][bd_id] = 0.0
                if not bd_id in reactions_data['thermo']['global'].keys():
                    reactions_data['thermo']['global'][bd_id] = 0.0
                try:
                    if thermo_ceil>=brsynth_dict[bd_id]['value']>=thermo_floor:
                        norm_thermo = (brsynth_dict[bd_id]['value']-thermo_floor)/(thermo_ceil-thermo_floor)
                    elif brsynth_dict[bd_id]['value']<thermo_floor:
                        norm_thermo = 0.0
                    elif brsynth_dict[bd_id]['value']>thermo_ceil:
                        norm_thermo = 1.0
                    norm_thermo = 1.0-norm_thermo
                    reactions_data['thermo']['reactions'][reac_id][bd_id] = norm_thermo
                    reactions_data['thermo']['global'][bd_id] += norm_thermo
                except (KeyError, TypeError) as e:
                    logging.warning('Cannot find the thermo: '+str(bd_id)+' for the reaction: '+str(reac_id))
        ####### FBA ##############
        #higher is better
        #return all the FBA values
        norm_fba = 0.0
        for bd_id in brsynth_dict:
            if bd_id[:4]=='fba_':
                reactions_data['fba']['reactions'][reac_id][bd_id] = 0.0
                try:
                    if fba_ceil>=brsynth_dict[bd_id]['value']>=fba_floor:
                        #min-max feature scaling
                        norm_fba = (brsynth_dict[bd_id]['value']-fba_floor)/(fba_ceil-fba_floor)
                    elif brsynth_dict[bd_id]['value']<=fba_floor:
                        norm_fba = 0.0
                    elif brsynth_dict[bd_id]['value']>fba_ceil:
                        norm_fba = 1.0
                    reactions_data['fba']['reactions'][reac_id][bd_id] = norm_fba
                except (KeyError, TypeError) as e:
                    logging.warning('Cannot find the objective: '+str(bd_id)+' for the reaction: '+str(reac_id))
    ##############################
    ##### target FBA value #######
    ##############################
    #higher is better
    #loop through all the different objectives and normalise the values
    #find the objective
    for obj_id in rpsbml_json['objectives']:
        brsynth_dict = rpsbml_json['objectives'][obj_id]
        if not obj_id in reactions_data['fba']['global'].keys():
            reactions_data['fba']['global'][obj_id] = 0.0
        try:
            if fba_ceil>=float(brsynth_dict['flux_value']['value'])>=fba_floor:
                norm_fba = (round(float(brsynth_dict['flux_value']['value']), 4)-fba_floor)/(fba_ceil-fba_floor)
            elif float(brsynth_dict['flux_value']['value'])<fba_floor:
                norm_fba = 0.0
            elif float(brsynth_dict['flux_value']['value'])>fba_ceil:
                norm_fba = 1.0
            reactions_data['fba']['global'][obj_id] = norm_fba
        except (KeyError, TypeError) as e:
            logging.warning('Could not retreive flux value: '+str(obj_id))
            continue
    try:
        target_norm_fba = reactions_data['fba']['global'][objective_id]
    except KeyError:
        target_norm_fba = 0.0
    brsynth_dict = rpsbml_json[pathway_id]
    for bd_id in brsynth_dict:
        if bd_id[:4]=='dfG_':
            try:
                reactions_data['thermo']['global'][bd_id] = reactions_data['thermo']['global'][bd_id]/float(members)
            except (KeyError, TypeError) as e:
                logging.warning('Cannot find the thermo: '+str(bd_id)+' for the reaction: '+str(bd_id))
    try:
        target_norm_thermo_mean = np.mean([reactions_data['thermo']['reactions'][i][thermo_id] for i in reactions_data['thermo']['reactions']])
    except KeyError:
        target_norm_thermo_mean = 0.0
    ############################
    ##### length of members ####
    ############################
    #lower is better
    if members>max_rp_steps:
        logging.warning('There are more steps than specified')
        norm_steps = 1.0
    else:
        try:
            norm_steps = (float(members)-1.0)/(float(max_rp_steps)-1.0)
        except ZeroDivisionError:
            norm_steps = 0.0
    norm_steps = 1.0-norm_steps
    ############################
    ####### rule_score #########
    ############################
    #higher is better
    #WARNING: using norm here is redundant but is used to detect the type of annotation
    norm_rScore = float(len(members))/reactions_data['rule_score']['global']
    rpsbml.addUpdateBRSynth(rp_pathway, 'norm_rule_scores', norm_rScore)
    ############################
    ##### global score #########
    ############################
    try:
        globalScore = (norm_steps*weight_rp_steps+
                       target_norm_fba*weight_fba+
                       norm_rScore*weight_rule_score+
                       target_norm_thermo_mean*weight_thermo
                       )/sum([weight_rp_steps, weight_fba, weight_thermo, weight_rule_score])
    except ZeroDivisionError:
        logging.warning('Global score calculation is dividing everything by 0')
        globalScore = 0.0
    rpsbml.addUpdateBRSynth(rp_pathway, 'global_score', globalScore)
    return globalScore
