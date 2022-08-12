Requires the "genomes" results from ATLAS (https://github.com/metagenome-atlas/atlas) this was written for ATLAS version 4.9.1

To use:
  1. Use a conda environment with ATLAS installed see https://github.com/metagenome-atlas/atlas for more information. 
  2. Append the "append_to_atlas_config_file" file to the atlas directory you wish to analyze with this package.
  3. You will have to put a symbolic link (or download the databases) to CAT/BAT taxonomy and database you would like to use see https://github.com/dutilh/CAT for details.
  5. You will have to put a symbolic link or the "CAT_Pack" folder from https://github.com/dutilh/CAT in the directory of the ATLAS results you will be analyzing. 
  6. while in the ATLAS directory you wish to analyze with this package run " {ATLAS conda environment} run all --snakefile {location of the Snakefile of this package} "
