import numpy as np
from openmmplumed import PlumedForce
from openmm.app import *
from openmm import *
from openmm.unit import *
import pdbfixer
import openmm_utils as op            # openmm functions
import sys 

def new_simulation(pdbf,nsteps,temp,pressure,use_plumed,plmdf,outfreq):
    positions,topology = op.fix_pdb(pdbf)                             # Fixes pdb with misiing atoms and terminal atoms

    #ff = ForceField('charmm36.xml', 'charmm36/water.xml')         # Creates forcefield
    ff = ForceField('amber99sbildn.xml', 'tip3p.xml')              # Creates forcefield - for DDR1 kinase

    positions,topology = op.add_hydrogen(positions,topology,ff)                          # Adds hydrogen according to the forcefield
    positions,topology = op.solvate_me(positions,topology,ff,True,1,'tip3p')     # (Default) neutralizes the charge with Na+ and Cl- ions

    #Create simulation
    simulation=op.get_LangevinM_system(topology,ff,temp,0.002)


    # Perform equilibration

    # Adding restrain
    simulation = op.add_pos_res(positions,topology,simulation,10)    # 10 KJ/(mol.A^2) contrain

    # Minimize
    positions,simulation = op.ener_Minimize(positions,simulation)

    # NVT - restrain
    positions,velocities,simulation=op.run_MD(positions,simulation,150000,ens='NVT',run_type='equilNVT')   #No '_' in run_type

    # NPT - restrain
    positions,velocities,simulation=op.run_MD(positions,simulation,150000,ens='NPT',run_type='equilNPT',\
                                           velocities=velocities,cont=True)                            #No '_' in run_type

    # Remove restrain
    simulation.context.setParameter('k',0.0*kilojoules_per_mole/angstroms**2)

#     # NVT - no restrain
#     positions,velocities,simulation=op.run_MD(positions,simulation,150000,ens='NVT',run_type='equilNVT',\
#                                            velocities=velocities,cont=True)                            #No '_' in run_type

    # NPT - no restrain
    positions,velocities,simulation=op.run_MD(positions,simulation,150000,ens='NPT',run_type='equilNPT',\
                                           velocities=velocities,cont=True)                            #No '_' in run_type


    # Run Simulations - NPT production run
    
    positions,velocities,simulation=op.run_MD(positions,simulation,nsteps,ens='NPT',run_type='prod',\
                                           writeXTC=True,XTC_information='protein',\
                                           use_plumed=use_plumed,plumed_file=plmdf,\
                                           velocities=velocities,cont=True,outfreq=outfreq)           #No '_' in run_type

def restart_simulation(pdbf,nsteps,temp,pressure,use_plumed,plmdf,outfreq):
    pdb=PDBFile(pdbf)
    positions=pdb.positions
    topology=pdb.topology
    #ff = ForceField('charmm36.xml', 'charmm36/water.xml')         # Creates forcefield
    ff = ForceField('amber99sbildn.xml', 'tip3p.xml')              # Creates forcefield - for DDR1 kinase

    simulation=op.get_LangevinM_system(topology,ff,temp,0.002)

    # Restart Simulations - NPT production run

    positions,velocities,simulation=op.run_MD(positions,simulation,nsteps,ens='NPT',run_type='prod_restart',\
                                           writeXTC=True,XTC_information='protein',\
                                           use_plumed=use_plumed,plumed_file=plmdf,\
                                           velocities=None,cont=False,restart=True,outfreq=outfreq)           #No '_' in run_type


if __name__=='__main__':
    '''
    Python file to perform Molecular dynamics using openmm - basic example to use openmm_utils package
    
        (A better script will be updated)
    
    Args:
    
    Input:
     -file            : PDB file to start 
     -nsteps          : Number simulation steps (default=25000000)
     -temp            : Temperature of simulation (default=300K)
     -pressure        : Pressure of simulation (default=1bar)
     -plumed          : Plumed file to bias or calculate CVs on the fly
    
    '''
    
    if '-file' in sys.argv:
        pdbf = sys.argv[sys.argv.index('-file')+1]
        flag=1
    else:
        print("An initial coordinate file is required")
        flag=0
    
    if '-nsteps' in sys.argv:
        nsteps = int(sys.argv[sys.argv.index('-nsteps')+1])
    else:
        nsteps = 25000000        # 100ns simulation                                
    
    if '-temp' in sys.argv:
        temp = int(sys.argv[sys.argv.index('-temp')+1])   #Integer value
    else:
        temp = 300              # 300K temperature
    
    if '-press' in sys.argv:
        pressure = float(sys.argv[sys.argv.index('-press')+1])   #float value
    else:
        pressure = 1.0              # 1.0 bar pressure
    
    #Should add an argument for pressure example as 1bar
        
    if '-plumed' in sys.argv:
        plmdf = sys.argv[sys.argv.index('-plumed')+1]
        use_plumed=True
    else:
        plmdf ='None'
        use_plumed=False
    
    if '-restart' in sys.argv:
        restart = sys.argv[sys.argv.index('-restart')+1]
    else:
        restart = False
    
    if '-outfreq' in sys.argv:
        outfreq = int(sys.argv[sys.argv.index('-outfreq')+1])   #Integer value
    else:
        outfreq = 1000  
        
        
        
    if flag:
        if not restart:
            new_simulation(pdbf,nsteps,temp,pressure,use_plumed,plmdf,outfreq)
            
        else:
            restart_simulation(pdbf,nsteps,temp,pressure,use_plumed,plmdf,outfreq)


