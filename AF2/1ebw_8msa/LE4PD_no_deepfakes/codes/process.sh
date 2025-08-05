#!/bin/bash -l
#SBATCH --job-name="process"
#SBATCH --output="process.log"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH --mem=50G
#SBATCH -t 48:00:00

#Note: this file has been setup to run on the shared nodes on the Comet supercomputer at the San Diego Supercomputing Center (SDSC).
#Documentation for this cluster is here:
#
#                                                 https://www.sdsc.edu/support/user_guides/comet.html
#
#Comet uses the SLURM scheduler to manage and submit jobs. To submit this job on Comet, run 'sbatch process.pbs'. 
#To run this file in the shell, use either 'sh process.sh' or 'bash process.sh'. If running locally, comment out the line
#reading 'module load gromacs'.
#
#This script also requires GROMACS (gromacs.org) be installed and parallelized with MPI. However, the processing commands here are
#run in serial, and thus a parallelized version of GROMACS is not required to run these processing commands. If the version of 
#GROMACS installed on your machine is not parallelized, simply change the line reading 'gmx=`which gmx_mpi`' to 'gmx=`which gmx`', which 
#will find the (serial) GROMACS executable on your machine.
#
#The topology file (the variable named 'top') can be any topology file recognized by GROMACS (.tpr, .pdb, .gro, etc.). Likewise with the trajectory ('traj')
#file. 
#
#Finally, this code is designed to write the alpha-carbons of the protein of interest to the .g96, .gro, and .pdb output files. The code assumes that the 
#group stored in index 1 of the topology file corresponds to the atoms of the protein and group 3 to the protein's alpha-carbons. If this is not the case,
#then the groups specified in the 'echo' commands prefixing the GROMACS processing commands must be adjusted accordingly -- the user should see the GROMACS
#documentation for more details. 

protname="pro" 
FILE_PATH="../prod/"
traj="${FILE_PATH}${protname}.xtc"
top="${FILE_PATH}${protname}.tpr"

if ! command -v gmx_mpi &> /dev/null
then 
	echo "gmx_mpi could not be found"
	
	if `command -v module load` &> /dev/null
	then
		module load gromacs
	fi

	if `command -v gmx_mpi` &> /dev/null
	then
		gmx=`command -v gmx_mpi`
	else
		gmx=`command -v gmx`
	fi
fi

echo $gmx
#Dump centered topology file
echo "1" | $gmx trjconv -f $traj -s $top -dump 0 -o top.pdb

$gmx editconf -f top.pdb -c -center 0 0 0 -o top.pdb

old_top=$top
top='top.pdb'
#Do this for PCA

# sasa
echo "1" | $gmx sasa -quiet -f $traj -s $top -or resarea.xvg -dt 10

#Correct PBCs
echo "3 1" | $gmx trjconv -quiet -f $traj -s $top -o after_pbc.xtc -pbc whole -center -boxcenter zero

echo "3 1" | $gmx trjconv -quiet -f after_pbc.xtc -s $top -fit rot -o after_rot.xtc

#Dump the processed trajectory, CA only
echo "3 3" | $gmx trjconv -quiet -f after_rot.xtc -s $top -o CA.xtc 
echo "3" | $gmx trjconv -quiet -f after_rot.xtc -s $top -o CA.pdb -dump 0 

exit

sh run_le4pd.sh './' $top 'CA.xtc'
