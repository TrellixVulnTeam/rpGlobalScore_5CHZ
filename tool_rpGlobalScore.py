#!/usr/bin/env python3

"""
Create on September 21 2019

@author: Melchior du Lac
@description: Galacxy specific script to call the function

python tool_rpGlobalScore.py -sbml test_rpThermo.rpsbml.xml -outputTar test_output.tar -weight_selenzyme 1.0 -weight_fba 1.0 -weight_thermo 1.0 -weight_rp_steps 1.0 -max_rp_steps 4 -topX 10 -thermo_ceil 8901.2 -thermo_floor -7570.2 -fba_ceil 999999.0 -fba_floor 0.0 -pathway_id rp_pathway 

"""

import sys
sys.path.insert(0, '/home/')

import rpToolServe
import argparse
import tempfile
import tarfile


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Given an SBML, extract the reaction rules and pass them to Selenzyme REST service and write the results to the SBML')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-sbml', type=str)
    parser.add_argument('-outputTar', type=str)
    parser.add_argument('-weight_selenzyme', type=float)
    parser.add_argument('-weight_fba', type=float)
    parser.add_argument('-weight_thermo', type=float)
    parser.add_argument('-weight_rp_steps', type=float)
    parser.add_argument('-max_rp_steps', type=int)
    parser.add_argument('-topX', type=int)
    parser.add_argument('-thermo_ceil', type=float)
    parser.add_argument('-thermo_floor', type=float)
    parser.add_argument('-fba_ceil', type=float)
    parser.add_argument('-fba_floor', type=float)
    parser.add_argument('-pathway_id', type=str)
    params = parser.parse_args()
    if params.sbml=='None' or params.sbml==None or params.sbml=='':
        if params.inputTar=='None' or params.inputTar==None or params.inputTar=='':
            logging.error('Cannot have no SBML and no TAR input')
            exit(0)
        rpToolServe.main(params.inputTar,
                         params.outputTar,
                         params.weight_rp_steps,
                         params.weight_selenzyme,
                         params.weight_fba,
                         params.weight_thermo,
                         params.max_rp_steps,
                         params.topX,
                         params.thermo_ceil,
                         params.thermo_floor,
                         params.fba_ceil,
                         params.fba_floor,
                         params.pathway_id)
    else:
        #make the tar.xz 
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            inputTar = tmpOutputFolder+'/tmp_input.tar.xz'
            with tarfile.open(inputTar, mode='w:xz') as tf:
                tf.add(params.sbml)
            rpToolServe.main(inputTar,
                             params.outputTar,
                             params.weight_rp_steps,
                             params.weight_selenzyme,
                             params.weight_fba,
                             params.weight_thermo,
                             params.max_rp_steps,
                             params.topX,
                             params.thermo_ceil,
                             params.thermo_floor,
                             params.fba_ceil,
                             params.fba_floor,
                             params.pathway_id)
