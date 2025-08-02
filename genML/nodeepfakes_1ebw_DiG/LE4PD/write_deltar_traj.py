#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import mdtraj as md
import sys

def load_traj(xtc, top, chunk_size = 5000, stride = 1):
	lst = []
	for chunk in md.iterload(xtc, chunk = chunk_size, stride = stride, top = top):
		print(chunk)
		lst.append(chunk)

	for i, dummy in enumerate(lst):
		print(dummy)
		NFRS = dummy.n_frames
		NATOMS = dummy.n_atoms
		if i == 0:
			traj =  np.reshape(dummy.xyz, (NFRS, 3*NATOMS))
		else:
			traj = np.vstack([traj, np.reshape(dummy.xyz, (NFRS, 3*NATOMS))])

	return traj.T, dummy


if len(sys.argv) < 3: raise Exception("Need two command line arguments: trajectory and topology files.")

traj = load_traj(sys.argv[1], sys.argv[2], stride = 1)
traj = traj[0]

#Get center-of-mass coordinates
for k in range(traj.shape[1]):
    traj[:,k] -= traj[:,k].mean(0)
      
dtraj = np.copy(traj)
dtraj = dtraj - (dtraj.mean(1)).reshape(dtraj.shape[0], 1)
for n, i in enumerate(range(0, dtraj.shape[0], 3)):
	np.savetxt('deltar_traj_' + str(n + 1), dtraj[i:i + 3,:])
    
