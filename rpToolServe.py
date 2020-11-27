#!/usr/bin/env python3

import libsbml
import argparse
import sys #exit using sys exit if any error is encountered
import os
import io
import logging
import tarfile
import glob
import json
import tempfile


sys.path.insert(0, '/home/')
import rpSBML
import rpTool

logging.basicConfig(
    #level=logging.DEBUG,
    #level=logging.WARNING,
    level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)


## run using HDD 3X less than the above function
#
#
#[0.10002239003499142, 0.13346271414277305, 0.6348436269211155, 0.13167126890112002]
def runGlobalScore_hdd(inputTar,
                       outputTar,
                       weight_rp_steps=0.10002239003499142,
                       weight_rule_score=0.13346271414277305,
                       weight_fba=0.6348436269211155,
                       weight_thermo=0.13167126890112002,
                       max_rp_steps=15,
                       topX=10,
                       thermo_ceil=5000.0,
                       thermo_floor=-5000.0,
                       fba_ceil=5.0,
                       fba_floor=0.0,
                       pathway_id='rp_pathway',
                       objective_id='obj_fraction',
                       thermo_id='dfG_prime_m'):
     """From a collection of rpSBML files, retreive the different characteristics of a pathway and combine them to calculate a global score.

    Note that the results are added to the passed dict directly.

    :param inputTar: An input tar collection of rpSBML files
    :param outputTar: An output tar collection of rpSBML files
    :param weight_rp_steps: The weight associated with the number of steps (Default: 0.10002239003499142)
    :param weight_rule_score: The weight associated with the mean of reaction rule scores (Default: 0.13346271414277305)
    :param weight_fba: The weight associated with the flux of the target (Default: 0.6348436269211155)
    :param weight_thermo: The weight associated with the sum of reaction Gibbs free energy (Default: 0.13167126890112002)
    :param max_rp_steps: The maximal number of steps are run in RP2 (Default: 15)
    :param thermo_ceil: The upper limit of Gibbs free energy for each reaction (Default: 5000.0)
    :param thermo_floor: The lower limit of Gibbs free energy for each reaction (Default: -5000.0)
    :param fba_ceil: The upper flux limit of the heterologous pathway (Default: 5.0)
    :param fba_floor: The lower flux limit of the heterologous pathway (Default: 5.0)
    :param pathway_id: The ID of the heterologous pathway (Default: rp_pathway)
    :param objective_id: The ID of the FBA objective (Default: obj_fraction)
    :param thermo_id: The ID of the Gibbs free energy that may be used. May be either dfG_prime_m or dfG_prime_o (Default: dfG_prime_m)

    :rtype: float
    :return: The global score
    """
    logging.debug('max_rp_steps: '+str(max_rp_steps))
    logging.debug('topX: '+str(topX))
    logging.debug('thermo_ceil: '+str(thermo_ceil))
    logging.debug('thermo_floor: '+str(thermo_floor))
    logging.debug('fba_ceil: '+str(fba_ceil))
    logging.debug('fba_floor: '+str(fba_floor))
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            if len(glob.glob(tmpInputFolder+'/*'))==0:
                logging.error('Input is empty')
                return {}
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
            top_file_names = [k for k, v in sorted(file_names_score.items(), key=lambda item: item[1], reverse=True)][:topX]
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpGlobalScore has not generated any results')
                return {}
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    file_name = str(sbml_path.split('/')[-1].replace('.rpsbml', '').replace('.sbml', '').replace('.xml', ''))
                    if file_name in top_file_names:
                        file_name += '.sbml.xml'
                        info = tarfile.TarInfo(file_name)
                        info.size = os.path.getsize(sbml_path)
                        ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return file_names_score
