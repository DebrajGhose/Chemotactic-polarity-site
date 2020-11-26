#!/bin/tcsh
#
#SBATCH -o screen.log
#SBATCH -e screen.err
#SBATCH  --mem=4G

module load Matlab/R2019a
matlab -nosplash -singleCompThread -r 'CopyFiles;quit'
