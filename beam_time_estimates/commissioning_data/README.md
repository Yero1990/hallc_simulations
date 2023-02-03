# Deuteron Commissioning Numerical Data

## Summary
This directory contains the numerical cross section data files for the Hall C Deuteron Commissioning Experiment (E12-10-003), for both data and various theoretical models which used the average kinematics input from the data to 
make the theoretical cross section calculations. 

### Directory Structure
`prl_data/` This directory contains the numerical cross sections (and reduced cross sections) for the overlappiong missing momenrum bins of the various datasets taken (80, 580\_set1, 580\_set2, 750\_set1, 750_set2, 750\_set3  ) for neutron recoil angle settings of 35, 45 and 75 +/- 45 deg.

`data/` This directory contains the final bin-centered, radiative-corrected numerical data cross sections as well as bin-centering correction factors, as well as the average kinematics corresponding to each bin. 

`theory/` This directory contains the raw theoretical numerical cross-sections  
from JVO (WJC2 potential), and MS (AV18 and CD-Bonn potentials). The JML (Paris potential) are actually in the `data/` directory, since these cross sections were calculated as part of the data analysis in SIMC, and were carried along with the experimental data numerical files.

`scripts/` This directory contains useful scripts for handling and plotting the above-mentioned numerical data files.


`numerical_data_combined/` This directory contains the numerical data directly extracted from our published PRL plots (fully combined values for overlapping bins), unlike the published numerical data, which have not been combined for overlapping bins.
