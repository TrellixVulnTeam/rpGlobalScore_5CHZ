#!/usr/bin/env python3

import libsbml
import sys
import logging

sys.path.insert(0, '/home/')
import rpSBML

## Extract the reaction SMILES from an SBML, query selenzyme and write the results back to the SBML
#
# NOTE: all the scores are normalised by their maximal and minimal, and normalised to be higher is better
#
def calculateGlobalScore(rpsbml,
                         weight_rp_steps,
                         weight_fba,
                         weight_thermo,
                         weight_reactionRule,
                         max_rp_steps,
                         pathway_id='rp_pathway',
                         rpFBAObj_name='rpFBA_obj'):
    groups = rpsbml.model.getPlugin('groups')
    rp_pathway = groups.getGroup(pathway_id)
    members = rp_pathway.getListOfMembers()
    all_rule_score = 0.0
    top_selenzyme = 0.0
    for member in members:
        reaction = rpsbml.model.getReaction(member.getIdRef())
        annot = reaction.getAnnotation()
        ####### rule_score ###########
        brsynth_dict = rpsbml.readBRSYNTHAnnotation(annot)
        all_rule_score += brsynth_dict['rule_score']
        ####### selenzyme ############
        try:
            #TODO: do a mean of the topX selenzyme score
            top_selenzyme += brsynth_dict['selenzyme'][sorted(brsynth_dict['selenzyme'], key=lambda kv: kv[1])[0]]
        except (IndexError, TypeError) as e:
            logging.warning('We are missing a selenzyme values')
    ##############################
    ####### rule_score ###########
    ##############################
    #lower is better
    norm_rule_score = all_rule_score/len(members)
    norm_rule_score = norm_rule_score-1.0
    #print('rule_score: '+str(norm_rule_score))
    ##############################
    ####### selenzyme ############
    ##############################
    #higher is better
    norm_selenzyme = (top_selenzyme/float(len(members)))/100.0
    #print('selenzyme: '+str(norm_selenzyme))
    ##############################
    ##### target FBA value ###
    ##############################
    #higher is better
    norm_fba = 0.0
    try:
        ###### theoretical flux limits
        #norm_fba = (brsynth_dict['fba_rpFBA_obj']['value'])/(999999.0)
        ###### practical flux limits
        #TODO: need to determine the best value of the maximal flux possible
        # look at the growth flux --> NO it dosn't work that way
        fbc = rpsbml.model.getPlugin('fbc')
        objective = fbc.getObjective(rpFBAObj_name)
        annot = objective.getAnnotation()
        brsynth_annot = annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
        for i in range(brsynth_annot.getNumChildren()):
            child = brsynth_annot.getChild(i)
            if child.getName()=='objective_value':
                #if the min is -999999.0 and max 999999.0 # ONLY use if reactions are reversible
                #norm_fba = (float(child.getAttrValue('value'))+999999.0)/(999999.0+999999.0)
                #if the min is 0.0 and max 999999.0
                norm_fba = round(float(child.getAttrValue('value')), 4)/999999.0
    except (KeyError, TypeError) as e:
        norm_fba = 0.0
    #print('FBA: '+str(norm_fba))
    ##############################
    #### group thermo ############
    ##############################
    #lower is better
    groups_list = groups.getListOfGroups()
    rp_group = groups_list[0]
    annot = rp_group.getAnnotation()
    brsynth_dict = rpsbml.readBRSYNTHAnnotation(annot)
    try:
        if 8901.2>brsynth_dict['dfG_prime_o']['value']>-7570.2:
            #if min is -7570.2 and max is 8901.2
            norm_thermo = (brsynth_dict['dfG_prime_o']['value']+7570.2)/(8901.2+7570.2)
            #if min is 8901.2 and max is -7570.2
            #WARNING: Does not work
            #norm_thermo = (brsynth_dict['dfG_prime_o']['value']+8901.2)/(-7570.2-8901.2)
        else:
            norm_thermo = 1.0
    except (KeyError, TypeError) as e:
        norm_thermo = 1.0
    norm_thermo = norm_thermo-1.0
    #print('Thermo: '+str(norm_thermo))
    ############################
    ##### length of members ####
    ############################
    #lower is better
    if len(members)>maxSetps:
        logging.warning('There are more steps than specified')
        norm_steps = 1.0
    else:
        norm_steps = (float(len(members))-1.0)/(float(max_rp_steps)-1.0)
    norm_steps = norm_steps-1.0
    #print('Steps: '+str(norm_steps))
    ############################
    ##### global score #########
    ############################
    #take the mean of the different normalised scores
    #globalScore = (norm_rule_score+norm_fba+norm_thermo+norm_selenzyme)/4.0
    globalScore = (norm_rule_score*weight_reactionRule+norm_steps*weight_rp_steps+norm_fba*weight_fba+norm_thermo*weight_thermo)/4.0
    #print('####### Global Score: '+str(globalScore)+' #########')
    #annot = groups.getAnnotation()
    bag_brsynth = annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
    #globalScore_child = bag_brsynth.getChild('globalScore')
    ##########################
    #### write annotation ####
    ##########################
    brsynth_annot = rp_pathway.getAnnotation().getChild('RDF').getChild('BRSynth').getChild('brsynth')
    tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:global_score value="'+str(globalScore)+'" /> </brsynth:brsynth>')
    brsynth_annot.addChild(tmpAnnot.getChild('global_score'))
