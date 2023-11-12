#!/bin/bash

analysis_script="analyze_simc_d2fsi.C"

#run_analysis="root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq28_fsi_rad\\\", 800, 28, \\\"pwia\\\" )\""
#echo $run_analysis

eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm500_thrq70_fsi_rad\\\",  500, 70, \\\"fsi\\\" )\""
eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm500_thrq70_pwia_rad\\\", 500, 70, \\\"pwia\\\" )\""


#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq28_fsi_rad\\\",  800, 28, \\\"fsi\\\" )\""
#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq28_pwia_rad\\\", 800, 28, \\\"pwia\\\" )\""

#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq49_fsi_rad\\\",  800, 49, \\\"fsi\\\" )\""
#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq49_pwia_rad\\\", 800, 49, \\\"pwia\\\" )\""

#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq55_fsi_rad\\\",  800, 55, \\\"fsi\\\" )\""
#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq55_pwia_rad\\\", 800, 55, \\\"pwia\\\" )\""

#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq60_fsi_rad\\\",  800, 60, \\\"fsi\\\" )\""
#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq60_pwia_rad\\\", 800, 60, \\\"pwia\\\" )\""

#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq66_fsi_rad\\\",  800, 66, \\\"fsi\\\" )\""
#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq66_pwia_rad\\\", 800, 66, \\\"pwia\\\" )\""

#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq72_fsi_rad\\\",  800, 72, \\\"fsi\\\" )\""
#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq72_pwia_rad\\\", 800, 72, \\\"pwia\\\" )\""

#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq79_fsi_rad\\\",  800, 79, \\\"fsi\\\" )\""
#eval "root -l -q -b  \"${analysis_script}( \\\"d2_pm800_thrq79_pwia_rad\\\", 800, 79, \\\"pwia\\\" )\""
