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


#sys.path.insert(0, '/home/')
import rpSBML
from rpTool import calculateGlobalScore


'''
##
#
#
def runGlobalScore_mem(inputTar, outputTar, weight_rp_steps, weight_fba, weight_thermo, weight_reactionRule, max_rp_steps, pathway_id, rpFBAObj_name):
    #loop through all of them and run FBA on them
    with tarfile.open(outputTar, 'w:xz') as tf:
        with tarfile.open(inputTar, 'r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(in_tf.extractfile(member).read().decode("utf-8")))
                    globalScore = calculateGlobalScore(rpsbml, weight_rp_steps, weight_fba, weight_thermo, weight_reactionRule, max_rp_steps, pathway_id, rpFBAObj_name)
                    data = libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
                    fiOut = io.BytesIO(data)
                    info = tarfile.TarInfo(member.name)
                    info.size = len(data)
                    tf.addfile(tarinfo=info, fileobj=fiOut)
'''


## run using HDD 3X less than the above function
#
#
def runGlobalScore_hdd(inputTar, outputTar, weight_rp_steps, weight_fba, weight_thermo, weight_reactionRule, max_rp_steps, topX, pathway_id, rpFBAObj_name):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(fileobj=inputTar, mode='r:xz')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            fileNames_score = {}
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '')
                rpsbml = rpSBML.rpSBML(fileName)
                rpsbml.readSBML(sbml_path)
                globalScore = calculateGlobalScore(rpsbml, weight_rp_steps, weight_fba, weight_thermo, weight_reactionRule, max_rp_steps, pathway_id, rpFBAObj_name)
                fileNames_score[fileName] = score
                rpsbml.writeSBML(tmpOutputFolder)
            #sort the results
            top_fileNames = [k for k, v in sorted(fileNames_score.items(), key=lambda item: item[1])][:topX]
            with tarfile.open(fileobj=outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', ''))
                    if fileName in top_fileNames:
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
                           float(params['weight_rp_steps']),
                           float(params['weight_fba']),
                           float(params['weight_thermo']),
                           float(params['weight_reactionRule']),
                           float(params['max_rp_steps']),
                           str(params['pathway_id']),
                           str(params['rpFBAObj_name']))
        '''
        #### HDD ####
        #weight_rp_steps, weight_fba, weight_thermo, pathway_id
        weight_reactionRule, max_rp_steps
        runGlobalScore_hdd(inputTar,
                           outputTar,
                           float(params['weight_rp_steps']),
                           float(params['weight_fba']),
                           float(params['weight_thermo']),
                           float(params['weight_reactionRule']),
                           float(params['max_rp_steps']),
                           int(params['topX']),
                           str(params['pathway_id']),
                           str(params['rpFBAObj_name']))
        ###### IMPORTANT ######
        outputTar.seek(0)
        #######################
        return send_file(outputTar, as_attachment=True, attachment_filename='rpGlobalScore.tar', mimetype='application/x-tar')


api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')


if __name__== "__main__":
    app.run(host="0.0.0.0", port=8888, debug=True, threaded=True)
