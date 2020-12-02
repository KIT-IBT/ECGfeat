#!/bin/bash

#SBATCH --ntasks=40
#SBATCH --time=0:30:00
#SBATCH --mem=64000

module load math/matlab
module load compiler/clang

matlab -batch evaluateFeatureAlgorithms
