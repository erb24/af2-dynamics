#!/bin/bash

which python
# $1 := PATH for locating files 
# $2 := topology file (PDB file)
# $3 := trajectory file (XTC file)
python run_le4pd.py $1 $2 $3
