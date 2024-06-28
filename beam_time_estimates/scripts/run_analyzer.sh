#!/bin/bash

# steering shell script to run multiple file analysis, read from an input file list
#analysis_script="analyze_simc_d2fsi.C"
#input_file='d2_fsi_basenames.txt'

analysis_script="analyze_simc_d2_test.C"
input_file='d2_pol_basenames.txt'

for line in $(cat $input_file)
do
    # check if line is commented out, if so then skip it
    if [ "${line:0:1}" == "#" ]; then
	continue
    else
	basename="$line"
	run_analysis="root -l -q -b  \"${analysis_script}( \\\"${basename}\\\")\""
	echo $run_analysis
	eval $run_analysis
	
    fi	
done
