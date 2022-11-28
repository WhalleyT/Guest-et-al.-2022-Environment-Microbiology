# Image analysis script

## Usage
This script takes csv files from the folder ```new_csvs``` and runs power estimates on them.

The outputs will be summarised in a csv file, and the values of each estimate with be stored in ```simulation_tmp```. ```environment.yml``` can be used to create a conda environment, ```job.sh``` is a template to run it on a slurm cluster.

To run from the command line (in Linux) ```Rscript power.R``` is sufficient.

## Citation
This code is in supplement to [Responses of Salmonella biofilms to oxidizing biocides: evidence of spatial clustering](https://sfamjournals.onlinelibrary.wiley.com/doi/full/10.1111/1462-2920.16263). [https://doi.org/10.1111/1462-2920.16263](https://doi.org/10.1111/1462-2920.16263).

## License

Shared under GPL. see ```License``` for more information.