#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --account=def-sievers
eval "$(conda shell.bash hook)"
conda activate pyplanet
python run.py
