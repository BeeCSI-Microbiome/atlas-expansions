Requires the "genomes" results from ATLAS (https://github.com/metagenome-atlas/atlas) this was written for ATLAS version 4.9.1

To use:
  1. Use a conda environment with ATLAS installed see https://github.com/metagenome-atlas/atlas for more information. 
  2. Append the "append_to_atlas_config_file" file to the atlas directory you wish to analyze with this package.  
  3. while in the ATLAS directory you wish to analyze with this package run " {ATLAS conda environment} run all --snakefile {location of the Snakefile of this package} "
