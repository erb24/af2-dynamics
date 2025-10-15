#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=4-24:15:40
#SBATCH --job-name=1ebw
#SBATCH --output=job_$SLURM_JOB_ID.out

source /data/aranganathana2/miniconda3/etc/profile.d/conda.sh
conda activate MD_py

python MD_openmm_adv_script2.py -file 1ebw.pdb -temp 300 -nsteps 500000000 -press 1 -outfreq 500

#python MD_openmm_adv.py -file NPT_prod_1.pdb -plumed WT_metad_restart_50.dat -temp 300 -nsteps 10000 -press 1
