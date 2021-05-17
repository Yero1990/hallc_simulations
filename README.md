# Hall C Simulations

This repository contains directories relevant to the
simulations (mainly using SIMC) done at Jefferson Lab Hall C.

### Directory Structure

* `utils/`: Relevant utility scripts mainly used to calculate the central kinematics of an experiment before simulations are done. 

* `simc_heep/`: Simulation data files/code used for hydrogen elastics *H(e,e'p)* reaction analysis

* `simc_deep/`: Simulation data files/code used for deuteron break-up *D(e,e'p)n* reaction analysis

* `infiles/`: Simulation input files required by SIMC (e.g. spectrometers, kinematics, acceptance, target, etc.)

* `outfiles/`: Simulation output files containing summary of the results. The three files have extension: .gen, .geni, .hist, and are generated for every simulation pass. Of particular importance in the normalization factor `normfac` in the .hist file which is necessary to properly account for the event weighting. See `deut_simc/simc.f` for definition of the normalization constant (tip: search for `normfac`)

* `worksim/`: Simulation output ASCII file containing the numerical quantities for all variables that will be written to a `ROOTTree`. The script fmake_tree.C on the parent directory (`../`) is used to convert the numerical quantities to a `ROOTFile` with the relevant leaf variables
in the TTree.

### Important Information 

If you plan to use this repository for doing hydrogen / deuteron simulations using SIMC, follow these instructions after cloning it for the initial setup:

1. Make sure the `deut_simc` (SIMC gfortran for deuteron analysis) has been cloned and compiled in `../hallc_simulations` by doing the following: 

 `git clone https://github.com/Yero1990/deut_simc`
    
2. Set up the binary file `grid.bin` read in by SIMC to do the proper deuteron cross-section weighting using the J.M. Laget theoretical cross sections (using Paris NN potential). You will need the .tar file which I can provide via e-mail if needed (email: cyero@jlab.org)

 * Put the .tar file in `../hallc_simulations` and do: `tar -zxvf Laget.zip`

  The `tar -zxvf` command will unpack the contents (i.e., cross-section numerical data files) onto a directory `Laget/deut_laget`  
  
  
3. Change back to the repository `hallc_simulations`  and execute the command:

  * `./set_symlinks.sh`

  to set up the proper symbolic links to the `deut_simc` repository as well as the  
  `deut_laget` directory so that the deuteron numerical cross-sections can be read during the
  simulation
    
### How-To Guide: Doing Simulations in SIMC
Once you have everything setup, then we can start doing simulations. This is a multi-step process which has been simplified by using a single steering script `simulate.py`. Simply include all the input files you want to simulate in `input_files.txt` and execute:

* `python simulate`

 
	
