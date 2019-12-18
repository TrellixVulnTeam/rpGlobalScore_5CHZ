#!/usr/bin/env python3

import requests
import argparse
import json


##
#
#
def rpGlobalScoreUpload(inputTar,
                        outputTar,
                        #weight_reactionRule,
                        weight_selenzyme,
                        weight_fba,
                        weight_thermo,
                        weight_rp_steps,
                        max_rp_steps,
                        topX,
                        pathway_id,
                        rpFBAObj_name,
                        server_url):
    # Post request
    data = {'pathway_id': pathway_id,
            'rpFBAObj_name': rpFBAObj_name,
            #'weight_reactionRule': weight_reactionRule,
            'weight_selenzyme': weight_selenzyme,
            'weight_fba': weight_fba,
            'weight_thermo': weight_thermo,
            'weight_rp_steps': weight_rp_steps,
            'max_rp_steps': float(max_rp_steps),
            'topX': topX}
    files = {'inputTar': open(inputTar, 'rb'),
             'data': ('data.json', json.dumps(data))}
    r = requests.post(server_url+'/Query', files=files)
    r.raise_for_status()
    with open(outputTar, 'wb') as ot:
        ot.write(r.content)


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
    parser.add_argument('-server_url', type=str)
    params = parser.parse_args()
    rpGlobalScoreUpload(params.inputTar,
                        params.outputTar,
                        #params.weight_reactionRule,
                        params.weight_selenzyme,
                        params.weight_fba,
                        params.weight_thermo,
                        params.weight_rp_steps,
                        params.max_rp_steps,
                        params.topX,
                        params.pathway_id,
                        params.rpFBAObj_name,
                        params.server_url)
