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
import glob
import shutil


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Given an SBML, extract the reaction rules and pass them to Selenzyme REST service and write the results to the SBML')
    parser.add_argument('-input', type=str)
    parser.add_argument('-input_format', type=str)
    parser.add_argument('-output', type=str)
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
    parser.add_argument('-obj_name', type=str)
    params = parser.parse_args()
    if params.input_format=='tar':
        if params.inputTar=='None' or params.inputTar==None or params.inputTar=='':
            logging.error('Cannot have no SBML and no TAR input')
            exit(0)
        rpToolServe.main(params.input,
                         params.output,
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
                         params.pathway_id,
                         params.obj_name)
    elif params.input_format=='sbml':
        #make the tar.xz 
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            input_tar = tmpOutputFolder+'/tmp_input.tar.xz'
            output_tar = tmpOutputFolder+'/tmp_output.tar.xz'
            with tarfile.open(input_tar, mode='w:xz') as tf:
                #tf.add(params.input)
                info = tarfile.TarInfo('single.rpsbml.xml') #need to change the name since galaxy creates .dat files
                info.size = os.path.getsize(params.input)
                tf.addfile(tarinfo=info, fileobj=open(params.input, 'rb')) 
            rpToolServe.main(input_tar,
                             output_tar,
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
                             params.pathway_id,
                             params.obj_name)
            with tarfile.open(output_tar) as outTar:
                outTar.extractall(tmpOutputFolder)
            out_file = glob.glob(tmpOutputFolder+'/*.rpsbml.xml')
            if len(out_file)>1:
                logging.warning('There are more than one output file...')
            shutil.copy(out_file[0], params.output)
    else:
        self.logging('Cannot identify the input_format: '+str(params.input_format))
