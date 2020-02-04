#!/usr/bin/env python3

import requests
import argparse
import json

"""
python tool_rpGlobalScore.py -input Galaxy832-[rpFBA].tar -output test_output.tar -weight_selenzyme 1.0 -weight_fba 1.0 -weight_thermo 1.0 -weight_rp_steps 1.0 -max_rp_steps 4 -topX 10 -thermo_ceil 8901.2 -thermo_floor -7570.2 -fba_ceil 999999.0 -fba_floor 0.0 -input_format tar -pathway_id rp_pathway -obj_name RP1_sink__restricted_biomass -server_url http://0.0.0.0:8881/REST

python tool_rpGlobalScore.py -input rp_1_1_rpFBA.rpsbml.xml -output test_output.rpsbml.xml -weight_selenzyme 1.0 -weight_fba 1.0 -weight_thermo 1.0 -weight_rp_steps 1.0 -max_rp_steps 4 -topX 10 -thermo_ceil 8901.2 -thermo_floor -7570.2 -fba_ceil 999999.0 -fba_floor 0.0 -input_format sbml -pathway_id rp_pathway -obj_name RP1_sink__restricted_biomass -server_url http://0.0.0.0:8881/REST

"""

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
                        objective_id,
                        thermo_id,
                        server_url):
    # Post request
    data = {'pathway_id': pathway_id,
            'obj_name': obj_name,	
            'weight_rp_steps': float(weight_rp_steps),
            'weight_selenzyme': float(weight_selenzyme),
            'weight_fba': float(weight_fba),
            'weight_thermo': float(weight_thermo),
            'max_rp_steps': float(max_rp_steps),
			'thermo_ceil': float(thermo_ceil),
			'thermo_floor': float(thermo_floor),
			'fba_ceil': float(fba_ceil),
			'fba_floor': float(fba_floor),
            'topX': int(topX),
            'objective_id': objective_id,
            'thermo_id': thermo_id}
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
    parser.add_argument('-objective_id', type=str)
    parser.add_argument('-thermo_id', type=str)
    parser.add_argument('-server_url', type=str)
    params = parser.parse_args()
    if params.input_format=='tar':
        rpGlobalScoreUpload(params.input,
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
                            params.objective_id,
                            params.thermo_id,
                            params.server_url)
    elif params.input_format=='sbml':
        #make the tar.xz 
        with tempfile.temporarydirectory() as tmpoutputfolder:
            input_tar = tmpoutputfolder+'/tmp_input.tar.xz'
            output_tar = tmpoutputfolder+'/tmp_output.tar.xz'
            with tarfile.open(input_tar, mode='w:xz') as tf:
                #tf.add(params.input)
                info = tarfile.tarinfo('single.rpsbml.xml') #need to change the name since galaxy creates .dat files
                info.size = os.path.getsize(params.input)
                tf.addfile(tarinfo=info, fileobj=open(params.input, 'rb')) 
        rpGlobalScoreUpload(input_tar,
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
                            params.objective_id,
                            params.thermo_id,
                            params.server_url)
            with tarfile.open(output_tar) as outtar:
                outtar.extractall(tmpoutputfolder)
            out_file = glob.glob(tmpoutputfolder+'/*.rpsbml.xml')
            if len(out_file)>1:
                logging.warning('there are more than one output file...')
            shutil.copy(out_file[0], params.output)
    else:
        self.logging('cannot identify the input_format: '+str(params.input_format))
