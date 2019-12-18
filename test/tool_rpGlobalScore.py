#!/usr/bin/env python3

"""
Create on September 21 2019

@author: Melchior du Lac
@description: Galacxy specific script to call the function

"""

import sys
sys.path.insert(0, '/home/')

import rpToolServe
import argparse


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Given an SBML, extract the reaction rules and pass them to Selenzyme REST service and write the results to the SBML')
    parser.add_argument('-inputTar', type=str)
    parser.add_argument('-outputTar', type=str)
    #parser.add_argument('-weight_reactionRule', type=float)
    parser.add_argument('-weight_selenzyme', type=float)
    parser.add_argument('-weight_fba', type=float)
    parser.add_argument('-weight_thermo', type=float)
    parser.add_argument('-weight_rp_steps', type=float)
    parser.add_argument('-max_rp_steps', type=int)
    parser.add_argument('-topX', type=int)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-rpFBAObj_name', type=str)
    params = parser.parse_args()
    rpToolServe.main(params.inputTar,
                     params.outputTar,
                     params.weight_rp_steps,
                     params.weight_selenzyme,
                     params.weight_fba,
                     params.weight_thermo,
                     #params.weight_reactionRule,
                     params.max_rp_steps,
                     params.topX,
                     params.pathway_id,
                     params.rpFBAObj_name)
