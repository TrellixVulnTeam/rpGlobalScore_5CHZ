#!/bin/bash

docker run -v ${PWD}/inside_run.sh:/home/inside_run.sh -v ${PWD}/tool_rpGlobalScore.py:/home/tool_rpGlobalScore.py -v ${PWD}/test_input.tar:/home/test_input.tar -v ${PWD}/results/:/home/results/ --rm brsynth/rpglobalscore /bin/sh /home/inside_run.sh

cp results/test_output.tar .
