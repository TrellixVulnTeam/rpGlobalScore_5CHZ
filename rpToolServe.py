#!/usr/bin/env python3

import libsbml
import argparse
import os
import io
import tarfile
import glob
import tempfile
import shutil

import rpSBML
import rpTool


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

##
#
#
def main(inputTar,
         outputTar,
         weight_rp_steps,
         weight_selenzyme,
         weight_fba,
         weight_thermo,
         max_rp_steps,
         topX,
         thermo_ceil=8901.2,
         thermo_floor=-7570.2,
         fba_ceil=999999.0,
         fba_floor=0.0,
         pathway_id='rp_pathway',
         objective_id='obj_rpFBA_frac',
         thermo_id='dfG_prime_m'):
    with open(inputTar, 'rb') as inputTar_bytes:
        outputTar_bytes = io.BytesIO()
        runGlobalScore_hdd(inputTar_bytes,
                           outputTar_bytes,
                           weight_rp_steps,
                           weight_selenzyme,
                           weight_fba,
                           weight_thermo,
                           max_rp_steps,
                           topX,
                           thermo_ceil,
                           thermo_floor,
                           fba_ceil,
                           fba_floor,
                           pathway_id,
                           objective_id,
                           thermo_id)
        ########## IMPORTANT #####
        outputTar_bytes.seek(0)
        ##########################
        with open(outputTar, 'wb') as f:
            shutil.copyfileobj(outputTar_bytes, f, length=131072)
