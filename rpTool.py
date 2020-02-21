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
                         weight_selenzyme,
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
                      'selenzyme': {'reactions': {}, 'global': 0.0},
                      'thermo': {'reactions': {}, 'global': {}}}
    #Loop through all the reactions
    for member in members:
        reactions_data['rule_score']['reactions'][member.getIdRef()] = 0.0
        reactions_data['fba']['reactions'][member.getIdRef()] = {}
        reactions_data['selenzyme']['reactions'][member.getIdRef()] = 0.0
        reactions_data['thermo']['reactions'][member.getIdRef()] = {}
        reaction = rpsbml.model.getReaction(member.getIdRef())
        annot = reaction.getAnnotation()
        brsynth_dict = rpsbml.readBRSYNTHAnnotation(annot)
        ####### selenzyme ############
        #higher is better
        try:
            #sum of mean of top selenzyme score
            top_selenzyme = brsynth_dict['selenzyme'][sorted(brsynth_dict['selenzyme'], key=lambda kv: kv[1])[0]]
            reactions_data['selenzyme']['reactions'][member.getIdRef()] = top_selenzyme/100.0
            reactions_data['selenzyme']['global'] += top_selenzyme/100.0
            rpsbml.addUpdateBRSynth(reaction, 'norm_selenzyme', top_selenzyme/100.0)
        except (IndexError, TypeError) as e:
            logging.warning('Missing selenzyme values for reaction: '+str(member.getIdRef()))
            rpsbml.addUpdateBRSynth(reaction, 'norm_selenzyme', 0.0)
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
                        #min-max feature scaling
                        norm_thermo = (brsynth_dict[bd_id]['value']-thermo_floor)/(thermo_ceil-thermo_floor)
                    elif brsynth_dict[bd_id]['value']<thermo_floor:
                        norm_thermo = 0.0
                    elif brsynth_dict[bd_id]['value']>thermo_ceil:
                        norm_thermo = 1.0
                    norm_thermo = 1.0-norm_thermo
                    reactions_data['thermo']['reactions'][member.getIdRef()][bd_id] = norm_thermo
                    reactions_data['thermo']['global'][bd_id] += norm_thermo
                    #if bd_id==thermo_id:
                    #    reactions_data['thermo']['global'][bd_id] += norm_thermo
                    rpsbml.addUpdateBRSynth(reaction, 'norm_'+bd_id, norm_thermo)
                except (KeyError, TypeError) as e:
                    logging.warning('Cannot find the thermo: '+str(bd_id)+' for the reaction: '+str(member.getIdRef()))
                    rpsbml.addUpdateBRSynth(reaction, 'norm_'+bd_id, 0.0)
                    #norm_thermo = 1.0
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
    #print('########### reactions_data ###########')
    #print(reactions_data)
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
        #print('------- objective ---------> '+str(objective.getId()))
        #print(brsynth_dict)
        try:
            #min-max feature scaling for FBA
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
        #update the flux obj
        for flux_obj in objective.getListOfFluxObjectives():
            brsynth_dict = rpsbml.readBRSYNTHAnnotation(flux_obj.getAnnotation())
            #min-max feature scaling
            try:
                norm_fba = (round(float(brsynth_dict['flux_value']['value']), 4)-fba_floor)/(fba_ceil-fba_floor)
            except (KeyError, TypeError) as e:
                norm_fba = 0.0
            rpsbml.addUpdateBRSynth(flux_obj, 'norm_flux_value', norm_fba)
    #print('@~@~@@~@~@~@~@~@~@~@')
    #print(reactions_data['fba']['global'])
    try:
        target_norm_fba = reactions_data['fba']['global'][objective_id]
    except KeyError:
        target_norm_fba = 0.0
    ##############################
    #### group thermo ############
    ##############################
    #lower is better
    #loop through all the reactions and retreive the thermo score for the reactions directly and sum up
    '''if you want to use the rpThermo pathway value
    brsynth_dict = rpsbml.readBRSYNTHAnnotation(rp_pathway.getAnnotation())
    norm_thermo = 0.0
    try:
        if thermo_ceil>=brsynth_dict['dfG_prime_m']['value']>=thermo_floor:
            #min-max feature scaling
            norm_thermo = (brsynth_dict['dfG_prime_m']['value']-thermo_floor)/(thermo_ceil-thermo_floor)
        elif brsynth_dict['dfG_prime_m']['value']<thermo_floor:
            norm_thermo = 0.0
        elif brsynth_dict['dfG_prime_m']['value']>thermo_ceil:
            norm_thermo = 1.0
    except (KeyError, TypeError) as e:
        norm_thermo = 1.0
    norm_thermo = 1.0-norm_thermo
    rpsbml.addUpdateBRSynth(rp_pathway, 'norm_dfG_prime_m', norm_thermo)
    '''
    '''
    #if you want the mean of each reaction normalisation
    for t_id in reactions_data['thermo']['global']:
        try:
            reactions_data['thermo']['global'][t_id] = reactions_data['thermo']['global'][t_id]/float(len(members))
            rpsbml.addUpdateBRSynth(rp_pathway, 'norm_dfG_prime_m', reactions_data['thermo']['global'][t_id])
        except KeyError:
            reactions_data['thermo']['global'][t_id] = 0.0
            rpsbml.addUpdateBRSynth(rp_pathway, 'norm_dfG_prime_m', reactions_data['thermo']['global'][t_id])
    '''
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
        #target_norm_thermo = reactions_data['thermo']['global'][thermo_id]
        target_norm_thermo_mean = np.mean([reactions_data['thermo']['reactions'][i][thermo_id] for i in reactions_data['thermo']['reactions']])
        target_norm_thermo_var = 1.0-stats.variation([reactions_data['thermo']['reactions'][i][thermo_id] for i in reactions_data['thermo']['reactions']])
    except KeyError:
        #target_norm_thermo = 0.0
        target_norm_thermo_mean = 0.0
        target_norm_thermo_var = 0.0
    ############################
    ##### length of members ####
    ############################
    #lower is better
    if len(members)>max_rp_steps:
        logging.warning('There are more steps than specified')
        norm_steps = 1.0
    else:
        try:
            norm_steps = (float(len(members))-1.0)/(float(max_rp_steps)-1.0)
        except ZeroDivisionError:
            norm_steps = 0.0
    norm_steps = 1.0-norm_steps
    rpsbml.addUpdateBRSynth(rp_pathway, 'norm_steps', norm_steps)
    ##############################
    ####### selenzyme ############
    ##############################
    #higher is better
    #NOTE: that missing selenzyme value will be considered 0.0 since we divide by the number of members
    reactions_data['selenzyme']['global'] = reactions_data['selenzyme']['global']/float(len(members))
    rpsbml.addUpdateBRSynth(rp_pathway, 'norm_selenzyme', reactions_data['selenzyme']['global'])
    #print('########### reactions_data global ###########')
    #print(reactions_data)
    ############################
    ##### global score #########
    ############################
    #globalScore = (norm_selenzyme*weight_selenzyme+norm_steps*weight_rp_steps+norm_fba*weight_fba+norm_thermo*weight_thermo)/4.0
    '''
    toDiv = 0.0
    if not weight_rp_steps==0.0:
        toDiv += 1.0
    if not weight_selenzyme==0.0:
        toDiv += 1.0
    if not weight_fba==0.0:
        toDiv += 1.0
    if not weight_thermo==0.0:
        toDiv += 1.0
    '''
    ##TODO: can check the thermo value and if + then linearly reduce its score
    #print('\t\t--------------------------- rpGlobalScore --------------------')
    #print('\t\tSelenzyme: '+str(reactions_data['selenzyme']['global'])+' * '+str(weight_selenzyme)+' = '+str(reactions_data['selenzyme']['global']*weight_selenzyme))
    #print('\t\tSteps: '+str(norm_steps)+' * '+str(weight_rp_steps)+' = '+str(norm_steps*weight_rp_steps))
    #print('\t\tFBA: '+str(target_norm_fba)+' * '+str(weight_fba)+' = '+str(target_norm_fba*weight_fba))
    #print('\t\tThermodynamics: '+str(target_norm_thermo)+' * '+str(weight_thermo)+' = '+str(target_norm_thermo*weight_thermo))
    try:
        globalScore = (reactions_data['selenzyme']['global']*weight_selenzyme+
                       norm_steps*weight_rp_steps+
                       target_norm_fba*weight_fba+
                       #target_norm_thermo*weight_thermo)/4.0
                       #target_norm_thermo*weight_thermo)/toDiv
                       (target_norm_thermo_mean*target_norm_thermo_var)*weight_thermo
                       )/sum([weight_rp_steps, weight_selenzyme, weight_fba, weight_thermo])
    except ZeroDivisionError:
        globalScore = 0.0
    #print('\t\tGlobal: '+str(globalScore))
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
                             weight_selenzyme,
                             weight_fba,
                             weight_thermo,
                             #weight_thermo_var,
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
                      'selenzyme': {'reactions': {}, 'global': 0.0},
                      'thermo': {'reactions': {}, 'global': {}}}
    members = len(rpsbml_json['reactions'])
    #Loop through all the reactions
    #reactions
    for reac_id in rpsbml_json['reactions']:
        reactions_data['rule_score']['reactions'][reac_id] = 0.0
        reactions_data['fba']['reactions'][reac_id] = {}
        reactions_data['selenzyme']['reactions'][reac_id] = 0.0
        reactions_data['thermo']['reactions'][reac_id] = {}
        brsynth_dict = rpsbml_json['reactions'][reac_id]
        ####### selenzyme ############
        #higher is better
        try:
            #sum of mean of top selenzyme score
            top_selenzyme = brsynth_dict['selenzyme'][sorted(brsynth_dict['selenzyme'], key=lambda kv: kv[1])[0]]
            reactions_data['selenzyme']['reactions'][reac_id] = top_selenzyme/100.0
            reactions_data['selenzyme']['global'] += top_selenzyme/100.0
        except (IndexError, TypeError) as e:
            logging.warning('Missing selenzyme values for reaction: '+str(reac_id))
        ####### Thermo ############
        #lower is better
        #WARNING: we will only take the dfG_prime_m value
        for bd_id in brsynth_dict:
            if bd_id[:4]=='dfG_':
                reactions_data['thermo']['reactions'][reac_id][bd_id] = 0.0
                if not bd_id in reactions_data['thermo']['global'].keys():
                    reactions_data['thermo']['global'][bd_id] = 0.0
                #reactions_data['thermo']['global'][bd_id] = 0.0
                try:
                    if thermo_ceil>=brsynth_dict[bd_id]['value']>=thermo_floor:
                        #min-max feature scaling
                        norm_thermo = (brsynth_dict[bd_id]['value']-thermo_floor)/(thermo_ceil-thermo_floor)
                    elif brsynth_dict[bd_id]['value']<thermo_floor:
                        norm_thermo = 0.0
                    elif brsynth_dict[bd_id]['value']>thermo_ceil:
                        norm_thermo = 1.0
                    norm_thermo = 1.0-norm_thermo
                    reactions_data['thermo']['reactions'][reac_id][bd_id] = norm_thermo
                    reactions_data['thermo']['global'][bd_id] += norm_thermo
                    #if bd_id==thermo_id:
                    #    reactions_data['thermo']['global'][bd_id] += norm_thermo
                except (KeyError, TypeError) as e:
                    logging.warning('Cannot find the thermo: '+str(bd_id)+' for the reaction: '+str(reac_id))
                    #norm_thermo = 1.0
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
    #print('########### reactions_data ###########')
    #print(reactions_data)
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
        #print('------- objective ---------> '+str(objective.getId()))
        #print(brsynth_dict)
        try:
            #min-max feature scaling for FBA
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
    #print('@~@~@@~@~@~@~@~@~@~@')
    #print(reactions_data['fba']['global'])
    try:
        target_norm_fba = reactions_data['fba']['global'][objective_id]
    except KeyError:
        target_norm_fba = 0.0
    ##############################
    #### group thermo ############
    ##############################
    #lower is better
    #loop through all the reactions and retreive the thermo score for the reactions directly and sum up
    '''if you want to use the rpThermo pathway value
    brsynth_dict = rpsbml.readBRSYNTHAnnotation(rp_pathway.getAnnotation())
    norm_thermo = 0.0
    try:
        if thermo_ceil>=brsynth_dict['dfG_prime_m']['value']>=thermo_floor:
            #min-max feature scaling
            norm_thermo = (brsynth_dict['dfG_prime_m']['value']-thermo_floor)/(thermo_ceil-thermo_floor)
        elif brsynth_dict['dfG_prime_m']['value']<thermo_floor:
            norm_thermo = 0.0
        elif brsynth_dict['dfG_prime_m']['value']>thermo_ceil:
            norm_thermo = 1.0
    except (KeyError, TypeError) as e:
        norm_thermo = 1.0
    norm_thermo = 1.0-norm_thermo
    rpsbml.addUpdateBRSynth(rp_pathway, 'norm_dfG_prime_m', norm_thermo)
    '''
    '''
    #if you want the mean of each reaction normalisation
    for t_id in reactions_data['thermo']['global']:
        try:
            reactions_data['thermo']['global'][t_id] = reactions_data['thermo']['global'][t_id]/float(len(members))
            rpsbml.addUpdateBRSynth(rp_pathway, 'norm_dfG_prime_m', reactions_data['thermo']['global'][t_id])
        except KeyError:
            reactions_data['thermo']['global'][t_id] = 0.0
            rpsbml.addUpdateBRSynth(rp_pathway, 'norm_dfG_prime_m', reactions_data['thermo']['global'][t_id])
    '''
    brsynth_dict = rpsbml_json[pathway_id]
    for bd_id in brsynth_dict:
        if bd_id[:4]=='dfG_':
            try:
                reactions_data['thermo']['global'][bd_id] = reactions_data['thermo']['global'][bd_id]/float(members)
            except (KeyError, TypeError) as e:
                logging.warning('Cannot find the thermo: '+str(bd_id)+' for the reaction: '+str(bd_id))
    try:
        #target_norm_thermo = reactions_data['thermo']['global'][thermo_id]
        target_norm_thermo_mean = np.mean([reactions_data['thermo']['reactions'][i][thermo_id] for i in reactions_data['thermo']['reactions']])
        target_norm_thermo_var = 1.0-stats.variation([reactions_data['thermo']['reactions'][i][thermo_id] for i in reactions_data['thermo']['reactions']])
    except KeyError:
        #target_norm_thermo = 0.0
        target_norm_thermo_mean = 0.0
        target_norm_thermo_var = 0.0
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
    ##############################
    ####### selenzyme ############
    ##############################
    #higher is better
    #NOTE: that missing selenzyme value will be considered 0.0 since we divide by the number of members
    reactions_data['selenzyme']['global'] = reactions_data['selenzyme']['global']/float(members)
    #print('########### reactions_data global ###########')
    #print(reactions_data)
    ############################
    ##### global score #########
    ############################
    #globalScore = (norm_selenzyme*weight_selenzyme+norm_steps*weight_rp_steps+norm_fba*weight_fba+norm_thermo*weight_thermo)/4.0
    '''
    toDiv = 0.0
    if not weight_rp_steps==0.0:
        toDiv += 1.0
    if not weight_selenzyme==0.0:
        toDiv += 1.0
    if not weight_fba==0.0:
        toDiv += 1.0
    if not weight_thermo==0.0:
        toDiv += 1.0
    '''
    ##TODO: can check the thermo value and if + then linearly reduce its score
    #print('\t\t--------------------------- rpGlobalScore --------------------')
    #print('\t\tSelenzyme: '+str(reactions_data['selenzyme']['global'])+' * '+str(weight_selenzyme)+' = '+str(reactions_data['selenzyme']['global']*weight_selenzyme))
    #print('\t\tSteps: '+str(norm_steps)+' * '+str(weight_rp_steps)+' = '+str(norm_steps*weight_rp_steps))
    #print('\t\tFBA: '+str(target_norm_fba)+' * '+str(weight_fba)+' = '+str(target_norm_fba*weight_fba))
    #print('\t\tThermodynamics: '+str(target_norm_thermo)+' * '+str(weight_thermo)+' = '+str(target_norm_thermo*weight_thermo))
    try:
        globalScore = (reactions_data['selenzyme']['global']*weight_selenzyme+
                       norm_steps*weight_rp_steps+
                       target_norm_fba*weight_fba+
                       #target_norm_thermo*weight_thermo)/4.0
                       #target_norm_thermo*weight_thermo)/toDiv
                       #target_norm_thermo*weight_thermo
                       #target_norm_thermo_mean*weight_thermo+
                       #target_norm_thermo_var*weight_thermo_var
                       (target_norm_thermo_mean*target_norm_thermo_var)*weight_thermo
                       )/sum([weight_rp_steps, weight_selenzyme, weight_fba, weight_thermo])
    except ZeroDivisionError:
        globalScore = 0.0
    return globalScore




