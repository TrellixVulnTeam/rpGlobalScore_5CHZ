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
                       weight_rp_steps=0.0,
                       weight_rule_score=0.5,
                       weight_fba=0.699707,
                       weight_thermo=0.8334961,
                       max_rp_steps=15,
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
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            file_names_score = {}
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                file_name = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
                rpsbml = rpSBML.rpSBML(file_name)
                rpsbml.readSBML(sbml_path)
                globalScore = rpTool.calculateGlobalScore_rpsbml(rpsbml,
                                                                 weight_rp_steps,
                                                                 weight_rule_score,
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
                file_names_score[file_name] = globalScore
                rpsbml.writeSBML(tmpOutputFolder)
            #sort the results
            top_file_names = [k for k, v in sorted(file_names_score.items(), key=lambda item: item[1])][:topX]
            with tarfile.open(outputTar, mode='w:xz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    file_name = str(sbml_path.split('/')[-1].replace('.rpsbml', '').replace('.sbml', '').replace('.xml', ''))
                    if file_name in top_file_names:
                        file_name += '.rpsbml.xml'
                        info = tarfile.TarInfo(file_name)
                        info.size = os.path.getsize(sbml_path)
                        ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return file_names_score

