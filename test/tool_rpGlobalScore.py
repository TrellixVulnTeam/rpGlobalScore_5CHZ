#!/usr/bin/env python3

import requests
import argparse
import json

##
#
#
def rpGlobalScoreUpload(inputTar,
                        outputTar,
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
                        obj_name,
                        server_url):
    # Post request
    data = {'pathway_id': pathway_id,
            'obj_name': obj_name,	
            'weight_rp_steps': weight_rp_steps,
            'weight_selenzyme': weight_selenzyme,
            'weight_fba': weight_fba,
            'weight_thermo': weight_thermo,
            'max_rp_steps': float(max_rp_steps),
			'thermo_ceil': float(thermo_ceil),
			'thermo_floor': float(thermo_floor),
			'fba_ceil': float(fba_ceil),
			'fba_floor': float(fba_floor),
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
    parser.add_argument('-weight_rp_steps', type=float)
    parser.add_argument('-weight_selenzyme', type=float)
    parser.add_argument('-weight_fba', type=float)
    parser.add_argument('-weight_thermo', type=float)
    parser.add_argument('-max_rp_steps', type=int)
    parser.add_argument('-topX', type=int)
    parser.add_argument('-thermo_ceil', type=float)
    parser.add_argument('-thermo_floor', type=float)
    parser.add_argument('-fba_ceil', type=float)
    parser.add_argument('-fba_floor', type=float)
    parser.add_argument('-pathway_id', type=str)
    parser.add_argument('-obj_name', type=str)
    parser.add_argument('-server_url', type=str)
    params = parser.parse_args()
    rpGlobalScoreUpload(params.inputTar,
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
                        params.pathway_id,
                        params.obj_name,
                        params.server_url)
