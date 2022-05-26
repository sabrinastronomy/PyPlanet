#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-sievers
eval "$(conda shell.bash hook)"
conda activate PyPlanet
python run.py
