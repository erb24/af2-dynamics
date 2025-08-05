#!/bin/bash

BD=$PWD
for seq in "3ttp" "2pc0" "1q9p" "1ebw" "4z4x" "6p9a"
do
	for method in 'BioEMU' 'AFc' 'DiG'
	do
		mkdir -v LE4PD
		cp -vr codes LE4PD
		cp -v resarea.xvg LE4PD
		cd LE4PD
		cp -v codes/run_le4pd.sh ./
		cp -v codes/*.py ./
		bash run_le4pd.sh ../ ../CA.pdb ../aligned.xtc
		cd $BD
	done
done
