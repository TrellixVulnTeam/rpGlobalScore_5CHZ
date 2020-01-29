#!/bin/bash

docker run -v ${PWD}/inside_run_sbml.sh:/home/inside_run_sbml.sh -v ${PWD}/tool_rpGlobalScore.py:/home/tool_rpGlobalScore.py -v ${PWD}/test_rpThermo.rpsbml.xml:/home/test_rpThermo.rpsbml.xml -v ${PWD}/results/:/home/results/ --rm brsynth/rpglobalscore /bin/sh /home/inside_run_sbml.sh

cp results/test_output.tar .
