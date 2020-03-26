#!/usr/bin/env python3

import libsbml
import argparse
import sys #exit using sys exit if any error is encountered
import os
import io
import tarfile
import glob
import json
import tempfile


sys.path.insert(0, '/home/')
import rpSBML
import rpTool


## run using HDD 3X less than the above function
#
#
def runGlobalScore_hdd(inputTar,
                       outputTar,
                       weight_rp_steps,
                       weight_fba,
                       weight_thermo,
                       max_rp_steps,
                       topX=10,
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
                        max_rp_steps,
                        topX=10,
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

