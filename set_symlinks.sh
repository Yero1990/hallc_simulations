#!/bin/bash

echo "Setting Up Symbolic Links . . ."

#local (my computer at home)
SIMC_DIR="../deut_simc/"


#Set up symbolic link for . . .

# 1) deuteron (h2) theory cross sections
ln -sf ${SIMC_DIR}h2.theory

# 2) parameter name lists (required to run simc from this directory)
ln -sf ${SIMC_DIR}nml_default.data

# 3) SIMC executable
ln -sf ${SIMC_DIR}simc

# 4) SIMC executable that takes input files from the command line
ln -sf ${SIMC_DIR}run_simc.sh

# 5) hms/  (necessary to read optics matrix)
ln -sf ${SIMC_DIR}hms

# 6) shms/ (necessary to read optics matrix)
ln -sf ${SIMC_DIR}shms

# 7) deut_laget for the numerical cross sections to be read in by SIMC
#  (Assuming the Laget.tar file has been unpacked on upper level directory (../) )
ln -sf "../Laget/deut_laget"

# Depending on the simulation output file size (put in ./worksim)
# this directory may have to be symbolically linked to a location where
# more space can be allocated. For example, if working on ifarm, /worksim
# may need to be linked to some other location is file size gets too big.
# On local computer, I think it should not be a problem

# ln -sf ./path/to/worksim/
