#!/bin/bash

python tool_rpGlobalScore.py -inputTar test_input.tar -outputTar test_output.tar -pathway_id rp_pathway -weight_selenzyme 1.0 -weight_fba 1.0 -weight_thermo 1.0 -weight_rp_steps 1.0 -max_rp_steps 15 -topX 10 -rpFBAObj_name rpFBA_obj
mv test_output.tar results/
