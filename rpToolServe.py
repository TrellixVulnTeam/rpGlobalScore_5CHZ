#!/usr/bin/env python3

import sys #exit using sys exit if any error is encountered
sys.path.insert(0, '/home/')

import libsbml
import argparse
import os
import io
import tarfile
import glob
import tempfile
import shutil

import rpSBML
from rpTool import calculateGlobalScore


'''
##
#
#
def runGlobalScore_mem(inputTar_bytes,
                       outputTar_bytes,
                       weight_rp_steps,
                       weight_selenzyme,
                       weight_fba,
                       weight_thermo,
                       weight_reactionRule,
                       max_rp_steps,
                       pathway_id,
                       rpFBAObj_name):
    #loop through all of them and run FBA on them
    with tarfile.open(outputTar_bytes, 'w:xz') as tf:
        with tarfile.open(inputTar_bytes, 'r:xz') as in_tf:
            for member in in_tf.getmembers():
                if not member.name=='':
                    rpsbml = rpSBML.rpSBML(member_name, libsbml.readSBMLFromString(in_tf.extractfile(member).read().decode("utf-8")))
                    globalScore = calculateGlobalScore(rpsbml, weight_rp_steps, weight_selenzyme, weight_fba, weight_thermo, weight_reactionRule, max_rp_steps, pathway_id, rpFBAObj_name)
                    data = libsbml.writeSBMLToString(rpsbml.document).encode('utf-8')
                    fiOut = io.BytesIO(data)
                    info = tarfile.TarInfo(member.name)
                    info.size = len(data)
                    tf.addfile(tarinfo=info, fileobj=fiOut)
'''


## run using HDD 3X less than the above function
#
#
def runGlobalScore_hdd(inputTar_bytes,
                       outputTar_bytes,
                       weight_rp_steps,
                       weight_selenzyme,
                       weight_fba,
                       weight_thermo,
                       #weight_reactionRule,
                       max_rp_steps,
                       topX,
                       pathway_id,
                       rpFBAObj_name):
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(fileobj=inputTar_bytes, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            fileNames_score = {}
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '')
                rpsbml = rpSBML.rpSBML(fileName)
                rpsbml.readSBML(sbml_path)
                globalScore = calculateGlobalScore(rpsbml,
                                                   weight_rp_steps,
                                                   weight_selenzyme,
                                                   weight_fba,
                                                   weight_thermo,
                                                   #weight_reactionRule,
                                                   max_rp_steps,
                                                   pathway_id,
                                                   rpFBAObj_name)
                fileNames_score[fileName] = globalScore
                rpsbml.writeSBML(tmpOutputFolder)
            #sort the results
            top_fileNames = [k for k, v in sorted(fileNames_score.items(), key=lambda item: item[1])][:topX]
            with tarfile.open(fileobj=outputTar_bytes, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', ''))
                    if fileName in top_fileNames:
                        fileName += '.sbml.xml'
                        info = tarfile.TarInfo(fileName)
                        info.size = os.path.getsize(sbml_path)
                        ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))

##
#
#
def main(inputTar,
         outputTar,
         weight_rp_steps,
         weight_selenzyme,
         weight_fba,
         weight_thermo,
         #weight_reactionRule,
         max_rp_steps,
         topX,
         pathway_id,
         rpFBAObj_name):
    with open(inputTar, 'rb') as inputTar_bytes:
        #print(inputTar_bytes.read())
        outputTar_bytes = io.BytesIO()
        runGlobalScore_hdd(inputTar_bytes,
                           outputTar_bytes,
                           weight_rp_steps,
                           weight_selenzyme,
                           weight_fba,
                           weight_thermo,
                           #weight_reactionRule,
                           max_rp_steps,
                           topX,
                           pathway_id,
                           rpFBAObj_name)
        ########## IMPORTANT #####
        outputTar_bytes.seek(0)
        ##########################
        with open(outputTar, 'wb') as f:
            shutil.copyfileobj(outputTar_bytes, f, length=131072)
