#!/usr/bin/env python3

import libsbml
import argparse
import sys #exit using sys exit if any error is encountered
import os
import io
import tarfile
import glob
import tempfile

import json
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort
from flask_restful import Resource, Api
import tempfile

sys.path.insert(0, '/home/')
import rpSBML

## Extract the reaction SMILES from an SBML, query selenzyme and write the results back to the SBML
#
#
def calculateGlobalScore(rpsbml, weight_steps, weight_fba, weight_thermo, pathway_id='rp_pathway', rpFBAObj_name='rpFBA_obj'):
    groups = rpsbml.model.getPlugin('groups')
    rp_pathway = groups.getGroup(pathway_id)
    members = rp_pathway.getListOfMembers()
    all_rule_score = 0.0
    all_selenzyme = 0.0
    for member in members:
        reaction = rpsbml.model.getReaction(member.getIdRef())
        annot = reaction.getAnnotation()
        ##############################
        ####### rule_score ###########
        ##############################
        brsynth_dict = rpsbml.readBRSYNTHAnnotation(annot)
        all_rule_score += brsynth_dict['rule_score']
        ##############################
        ####### selenzyme ############
        ##############################
        try:
            #top selenzyme score
            all_selenzyme += brsynth_dict['selenzyme'][sorted(brsynth_dict['selenzyme'], key=lambda kv: kv[1])[0]]
        except (IndexError, TypeError) as e:
            logging.warning('We are missing a selenzyme values')
    #higher is better
    norm_selenzyme = (all_selenzyme/float(len(members)))/100.0
    #print('selenzyme: '+str(norm_selenzyme))
    #lower is better
    norm_rule_score = all_rule_score/len(members)
    norm_rule_score = 1.0-norm_rule_score
    #print('rule_score: '+str(norm_rule_score))
    ##############################
    ##### target FBA value ###
    ##############################
    #higher is better
    #fbc_obj = fbc_p.getObjective(rpFBAObj_name)
    #reaction = rpsbml.model.getReaction(obj.getListOfFluxObjectives()[0].getReaction())
    #annot = reaction.getAnnotation()
    #brsynth_dict = rpsbml.readBRSYNTHAnnotation(annot)
    norm_fba = 0.0
    try:
        ###### theoretical flux limits
        #if the min is -999999.0 and max 999999.0
        #norm_fba = (brsynth_dict['fba_rpFBA_obj']['value']+999999.0)/(999999.0+999999.0)
        #if the min is 0.0 and max 999999.0
        #norm_fba = (brsynth_dict['fba_rpFBA_obj']['value'])/(999999.0)
        ###### practical flux limits
        #TODO: need to determine the best value of the maximal flux possible
        # look at the growth flux
        fbc = rpsbml.model.getPlugin('fbc')
        objective = fbc.getObjective('rpFBA_obj')
        annot = objective.getAnnotation()
        brsynth_annot = annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
        for i in range(brsynth_annot.getNumChildren()):
            child = brsynth_annot.getChild(i)
            if child.getName()=='objective_value':
                norm_fba = (float(child.getAttrValue('value')))/(100.0)
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
        #if min is -7570.2 and max is 8901.2
        #norm_thermo = (brsynth_dict['dfG_prime_o']['value']-8901.2)/(-7570.2+8901.2)
        #if min is 8901.2 and max is -7570.2
        norm_thermo = (brsynth_dict['dfG_prime_o']['value']+7570.2)/(8901.2+7570.2)
    except (KeyError, TypeError) as e:
        norm_thermo = (7570.2)/(8901.2+7570.2)
    norm_thermo = 1.0-norm_thermo
    #print('Thermo: '+str(norm_thermo))
    ############################
    ##### length of members ####
    ############################
    #lower is better
    norm_steps = (float(len(members))-1.0)/(15.0-1.0)
    norm_steps = 1.0-norm_steps
    #print('Steps: '+str(norm_steps))
    ############################
    ##### global score #########
    ############################
    #take the mean of the different normalised scores
    #globalScore = (norm_rule_score+norm_fba+norm_thermo+norm_selenzyme)/4.0
    globalScore = (norm_steps*weight_steps+norm_fba*weight_fba+norm_thermo*weight_thermo)/3.0
    #print('####### Global Score: '+str(globalScore)+' #########')
    #annot = groups.getAnnotation()
    bag_brsynth = annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
    #globalScore_child = bag_brsynth.getChild('globalScore')
    #### annotation ####
    brsynth_annot = rp_pathway.getAnnotation().getChild('RDF').getChild('BRSynth').getChild('brsynth')
    tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<brsynth:brsynth xmlns:brsynth="http://brsynth.eu"> <brsynth:global_score value="'+str(globalScore)+'" /> </brsynth:brsynth>')
    brsynth_annot.addChild(tmpAnnot.getChild('global_score'))



##
#
#
def runGlobalScore_mem(inputTar, outputTar, weight_steps, weight_fba, weight_thermo, pathway_id):
    #loop through all of them and run FBA on them
    with tarfile.open(outputTar, 'w:xz') as tf:
        with tarfile.open(inputTar, 'r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(in_tf.extractfile(member).read().decode("utf-8")))
                    calculateGlobalScore(rpsbml, weight_steps, weight_fba, weight_thermo, pathway_id)
                    data = libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
                    fiOut = io.BytesIO(data)
                    info = tarfile.TarInfo(member.name)
                    info.size = len(data)
                    tf.addfile(tarinfo=info, fileobj=fiOut)


## run using HDD 3X less than the above function
#
#
def runGlobalScore_hdd(inputTar, outputTar, weight_steps, weight_fba, weight_thermo, pathway_id):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(inputTar, 'r:xz')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '')
                rpsbml = rpSBML.rpSBML(fileName)
                rpsbml.readSBML(sbml_path)
                calculateGlobalScore(rpsbml, weight_steps, weight_fba, weight_thermo, pathway_id)
                rpsbml.writeSBML(tmpOutputFolder)
            with tarfile.open(outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', ''))
                    fileName += '.sbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))


#######################################################
############## REST ###################################
#######################################################


app = Flask(__name__)
api = Api(app)


def stamp(data, status=1):
    appinfo = {'app': 'rpGlobalScore', 'version': '1.0',
               'author': 'Melchior du Lac',
               'organization': 'BRS',
               'time': datetime.now().isoformat(),
               'status': status}
    out = appinfo.copy()
    out['data'] = data
    return out


class RestApp(Resource):
    """ REST App."""
    def post(self):
        return jsonify(stamp(None))
    def get(self):
        return jsonify(stamp(None))


class RestQuery(Resource):
    """ REST interface that generates the Design.
        Avoid returning numpy or pandas object in
        order to keep the client lighter.
    """
    def post(self):
        inputTar = request.files['inputTar']
        params = json.load(request.files['data'])
        #pass the files to the rpReader
        outputTar = io.BytesIO()
        #### MEM ####
        '''
        runGlobalScore_mem(inputTar,
                           outputTar,
                           int(params['weight_steps']),
                           int(params['weight_fba']),
                           int(params['weight_thermo']))
        '''
        #### HDD ####
        #weight_steps, weight_fba, weight_thermo, pathway_id
        runGlobalScore_hdd(inputTar,
                           outputTar,
                           int(params['weight_steps']),
                           int(params['weight_fba']),
                           int(params['weight_thermo']))
        ###### IMPORTANT ######
        outputTar.seek(0)
        #######################
        return send_file(outputTar, as_attachment=True, attachment_filename='rpGlobalScore.tar', mimetype='application/x-tar')


api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')


if __name__== "__main__":
    app.run(host="0.0.0.0", port=8888, debug=True, threaded=True)
