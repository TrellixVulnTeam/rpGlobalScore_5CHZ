#!/bin/sh

docker run -d -p 8888:8888 --name test_rpGlobalScore brsynth/rpglobalscore
sleep 10
python tool_rpGlobalScore.py -inputTar test_input.tar -outputTar test_output.tar -weight_rp_steps 1.0 -weight_selenzyme 1.0 -weight_fba 1.0 -weight_thermo 1.0 -max_rp_steps 15 -topX 10 -pathway_id rp_pathway -rpFBAObj_name rpFBA_obj -server_url http://0.0.0.0:8888/REST
docker kill test_rpGlobalScore
docker rm test_rpGlobalScore
