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
import rpTool

'''
##
#
#
def runGlobalScore_mem(inputTar,
                       outputTar,
                       weight_rp_steps,
                       weight_selenzyme,
                       weight_fba,
                       weight_thermo,
                       #weight_reactionRule,
                       max_rp_steps,
                       pathway_id,
                       obj_name):
    #loop through all of them and run FBA on them
    with tarfile.open(outputTar, 'w:xz') as tf:
        with tarfile.open(inputTar, 'r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(in_tf.extractfile(member).read().decode("utf-8")))
                    globalScore = calculateGlobalScore(rpsbml,
                                                       weight_rp_steps,
                                                       weight_selenzyme,
                                                       weight_fba,
                                                       weight_thermo,
                                                       #weight_reactionRule,
                                                       max_rp_steps,
                                                       pathway_id,
                                                       obj_name)
                    data = libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
                    fiOut = io.BytesIO(data)
                    info = tarfile.TarInfo(member.name)
                    info.size = len(data)
                    tf.addfile(tarinfo=info, fileobj=fiOut)
'''


## run using HDD 3X less than the above function
#
#
def runGlobalScore_hdd(inputTar,
                       outputTar,
                       weight_rp_steps,
                       weight_selenzyme,
                       weight_fba,
                       weight_thermo,
                       max_rp_steps,
                       topX,
                       thermo_ceil=8901.2,
                       thermo_floor=-7570.2,
                       fba_ceil=5.0,
                       fba_floor=0.0,
                       pathway_id='rp_pathway',
                       objective_id='obj_RP1_sink__restricted_biomass',
                       thermo_id='dfG_prime_m'):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(fileobj=inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            fileNames_score = {}
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                rpsbml = rpSBML.rpSBML(fileName)
                rpsbml.readSBML(sbml_path)
                globalScore = rpTool.calculateGlobalScore(rpsbml,
                                                   weight_rp_steps,
                                                   weight_selenzyme,
                                                   weight_fba,
                                                   weight_thermo,
                                                   max_rp_steps,
                                                   thermo_ceil,
                                                   thermo_floor,
                                                   fba_ceil,
                                                   fba_floor,
                                                   pathway_id,
                                                   objective_id,
                                                   thermo_id)
                fileNames_score[fileName] = globalScore
                rpsbml.writeSBML(tmpOutputFolder)
            #sort the results
            top_fileNames = [k for k, v in sorted(fileNames_score.items(), key=lambda item: item[1])][:topX]
            with tarfile.open(fileobj=outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.rpsbml', '').replace('.sbml', '').replace('.xml', ''))
                    if fileName in top_fileNames:
                        fileName += '.rpsbml.xml'
                        info = tarfile.TarInfo(fileName)
                        info.size = os.path.getsize(sbml_path)
                        ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return fileNames_score

## RetroPath2.0 reader for local packages
#
#
def runGlobalScore_json(rpsbml_json,
                       weight_rp_steps,
                       weight_selenzyme,
                       weight_fba,
                       weight_thermo,
                       #weight_thermo_var,
                       max_rp_steps,
                       topX,
                       thermo_ceil=8901.2,
                       thermo_floor=-7570.2,
                       fba_ceil=5.0,
                       fba_floor=0.0,
                       pathway_id='rp_pathway',
                       objective_id='obj_RP1_sink__restricted_biomass',
                       thermo_id='dfG_prime_m'):
    fileNames_score = {}
    for rpsbml_id in rpsbml_json:
        globalScore = rpTool.calculateGlobalScore_json(rpsbml_json[rpsbml_id],
                                           weight_rp_steps,
                                           weight_selenzyme,
                                           weight_fba,
                                           weight_thermo,
                                           #weight_thermo_var,
                                           max_rp_steps,
                                           thermo_ceil,
                                           thermo_floor,
                                           fba_ceil,
                                           fba_floor,
                                           pathway_id,
                                           objective_id,
                                           thermo_id)
        fileNames_score[rpsbml_id] = globalScore
    return fileNames_score

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
                           float(params['weight_rp_steps']),
                           float(params['weight_selenzyme']),
                           float(params['weight_fba']),
                           float(params['weight_thermo']),
                           #float(params['weight_reactionRule']),
                           float(params['max_rp_steps']),
                           str(params['pathway_id']),
                           str(params['obj_name']))
        '''
        #### HDD ####
        #weight_rp_steps, weight_fba, weight_thermo, pathway_id
        runGlobalScore_hdd(inputTar,
                           outputTar,
                           float(params['weight_rp_steps']),
                           float(params['weight_selenzyme']),
                           float(params['weight_fba']),
                           float(params['weight_thermo']),
                           float(params['max_rp_steps']),
                           int(params['topX']),
                           float(params['thermo_ceil']),
                           float(params['thermo_floor']),
                           float(params['fba_ceil']),
                           float(params['fba_floor']),
                           str(params['pathway_id']),
                           str(params['objective_id']),
                           str(params['thermo_id']))
        ###### IMPORTANT ######
        outputTar.seek(0)
        #######################
        return send_file(outputTar, as_attachment=True, attachment_filename='rpGlobalScore.tar', mimetype='application/x-tar')


api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')


if __name__== "__main__":
    app.run(host="0.0.0.0", port=8888, debug=False, threaded=True)
