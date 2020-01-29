#!/bin/bash

python tool_rpGlobalScore.py -sbml test_rpThermo.rpsbml.xml -outputTar test_output.tar -weight_selenzyme 1.0 -weight_fba 1.0 -weight_thermo 1.0 -weight_rp_steps 1.0 -max_rp_steps 4 -topX 10 -thermo_ceil 8901.2 -thermo_floor -7570.2 -fba_ceil 999999.0 -fba_floor 0.0 -pathway_id rp_pathway 
mv test_output.tar results/
