#!/bin/bash 
for i in `seq 1 71`
do
  cp -v tcf_DUMMY.pbs tcf_${i}.pbs
  sed -i "s#DUMMY#${i}#g" tcf_${i}.pbs
  cp -v tcf_DUMMY.f95 ./tcf_${i}.f95
  sbatch tcf_${i}.pbs
done

exit
