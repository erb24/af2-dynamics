from openmm import *
from openmm.app import *
from openmm.unit import *
import mdtraj as md
from mdtraj.reporters import XTCReporter 
import re
import pdbfixer
from openmmplumed import PlumedForce
import time
import numpy as np

################################# File Management ###################################

def check_file(fname):
    '''
    Checks the existence of a file in the given pathway and gives a newfile name
    filename format  - <string_indentifier>_<int_identifier>.<extension>
    example filename - fixed_1.pdb
    '''
    if os.path.isfile(fname):
        ident=fname.split('.')[0].split('_')
        fname=f"{'_'.join([wrd for wrd in ident[0:-1]])}_{int(ident[-1])+1}.{fname.split('.')[1]}"
        fname=check_file(fname)
    return fname

############################ Structure Editing Functions ############################

# Will be added into a class soon!

#structuture
def fix_pdb(file_n):
    """
    fixes the raw pdb from colabfold using pdbfixer.
    This needs to be performed to cleanup the pdb and to start simulation 
    Fixes performed: missing residues, missing atoms and missing Terminals
    """
    raw_pdb=file_n;

    # fixer instance
    fixer = pdbfixer.PDBFixer(raw_pdb)

    #finding and adding missing residues including terminals
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    outfile=check_file('fixed_1.pdb')
    PDBFile.writeFile(fixer.topology, fixer.positions, open(outfile,'w'), keepIds=True)
    return fixer.positions, fixer.topology


def add_hydrogen(positions,topology,forcefield,write_file=True,):
    """
    Adds missing hydrogen to the pdb for a particular forcefield
    """
    modeller = Modeller(topology, positions)
    modeller.addHydrogens(forcefield);
    if write_file:
        hydfile=check_file('fixedH_1.pdb')
        PDBFile.writeFile(modeller.topology, modeller.positions, open(hydfile, 'w'))
    return modeller.positions, modeller.topology

def solvate_me(positions,topology,forcefield,write_file=True,padding=1,\
               water_model='tip3p',positiveIon='Na+',negativeIion='Cl-'):
    '''
    Creates a box of solvent with 1 nm padding and neutral charge
    '''
    modeller = Modeller(topology, positions)
    modeller.addSolvent(forcefield, padding=padding*nanometers, model=water_model, neutralize=True, positiveIon=positiveIon, negativeIon=negativeIion)
    if write_file:
        solvfile=check_file('solvated_1.pdb')
        PDBFile.writeFile(modeller.topology, modeller.positions, open(solvfile, 'w'))
    return modeller.positions, modeller.topology


############################ Perform Dynamics using openMM ############################

# Will be added into a class soon!

def get_LangevinM_system(topology,forcefield,temp=300,dt=0.002):
    system = forcefield.createSystem(topology,nonbondedMethod=PME,nonbondedCutoff=1.2*nanometer,switchDistance=1.0*nanometer,\
                                     constraints=HBonds);
    integrator = LangevinMiddleIntegrator(temp*kelvin, 1/picoseconds,dt*picoseconds);
    platform = Platform.getPlatformByName('CUDA');
    properties = {'Precision': 'double'} #change if required after setting up openmm
    simulation = Simulation(topology, system, integrator, platform)
    
    return simulation


def get_NoseHoover_system(topology,forcefield,temp=300,dt=0.002):
    system = forcefield.createSystem(topology,nonbondedMethod=PME,nonbondedCutoff=1.2*nanometer,switchDistance=1.0*nanometer,\
                                     constraints=HBonds);
    integrator = NoseHooverIntegrator(temp*kelvin, 1/picoseconds,dt*picoseconds);
    platform = Platform.getPlatformByName('CUDA');
    properties = {'Precision': 'double'} #change if required after setting up openmm
    simulation = Simulation(topology, system, integrator, platform)
    
    return simulation


def add_pos_res(positions,topology,simulation,k=1000):
    '''
    
    Adds an harmonic potential to the heavy atoms of the system(proteins) with an user defined force constant 'k'
    
    '''
    
    AA=['ALA','ASP','CYS','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ARG','PRO','GLN','ASN','SER','THR','VAL','TRP','TYR']
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2") # Harmonic potential for position restrain
    force.addGlobalParameter("k",k*kilocalories_per_mole/angstroms**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    
    index=0;
    for i, res in enumerate(topology.residues()):
        if res.name in AA:                              # Required to select only the protein atoms
            for at in res.atoms():
                if not re.search(r'H',at.name):         # All heavy Atoms -(exculdes Hydrogens)
                    force.addParticle(index,positions[index].value_in_unit(nanometers))
                index+=1;
                
    posres_sys=simulation.context.getSystem()                  # A gets System for a simulation instance
    posres_sys.addForce(force)                          # Modifies system with custom Force
    simulation.context.reinitialize()                          # initializes the simulation instance with the modified system
    
    return simulation


def ener_Minimize(positions,simulation,tolerance=10,n_iter=1500,write_file=True):
    '''
    Energy minimization step for simulation object
    '''
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy(tolerance=tolerance*kilojoule/mole,maxIterations=n_iter)
    minim_positions = simulation.context.getState(getPositions=True).getPositions()
    if write_file:
        minimfile=check_file('minim_1.pdb')
        pdb_positions = simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()
        PDBFile.writeFile(simulation.topology, pdb_positions, open(minimfile, 'w'))
    return minim_positions,simulation


def run_MD(positions,simulation,nsteps,ens='NVT',run_type='md',temp=300,pressure=1,use_plumed=False,plumed_file='None',\
           save_chkpt_file=True,outfreq=5000,writeXTC=True,XTC_information='protein',velocities=None,cont=False,restart=False):
    '''
    Function to run simulation for a particular ensemble. The function can also incorporate plumed file. 
    There is a restart option to start a simulation from checkpoint file
    There is also a continue option if multiple simulation is required to be run in a single script
    
    positions, simulation and nsteps are required parameters
    
    '''
    
    ## Add forces according to the ensemble or plumed file
    
    if ens=='NPT':
        simulation.context.getSystem().addForce(MonteCarloBarostat(pressure*bar, temp*kelvin)) # Pressure coupling
    
    if use_plumed:                                       # Adding custom force via plumed or getting CVs via plumed
        fid=open(plumed_file,'r')
        ff=fid.read()
        force=PlumedForce(ff)
        pl_system=simulation.context.getSystem()
        pl_system.addForce(force)
        simulation.context.reinitialize(True)         
    
    ## Initialize the simulation with positions and velocity (if you are continuing to run a simulation )
    
    simulation.context.setPositions(positions)
    if cont:
        simulation.context.setVelocities(velocities)
    
    
    ## If restarting from checkpoint file
    
    if restart:
        simulation.loadCheckpoint('chkptfile.chk')
    
    
    ## Append reporters for the simulation output and output files

    outfile=check_file(f'{ens}_{run_type}_1.pdb')
    logfile=check_file(f'{ens}_{run_type}_1.txt')          #output files
    xmlfile=check_file(f'{ens}_{run_type}_1.state')
    
    simulation.reporters=[]
    outlog=open(logfile,'w')
    simulation.reporters.append(StateDataReporter(outlog, outfreq*2, step=True,potentialEnergy=True,kineticEnergy=True,separator='\t|\t',progress=True,speed=True,totalSteps=nsteps))
    
    if save_chkpt_file:
        chkpt_freq=0.05*nsteps
        simulation.reporters.append(CheckpointReporter('chkptfile.chk', chkpt_freq))
    
    if writeXTC:                                     #default is True
        outfname=check_file(f'{ens}_{run_type}_1.xtc')
        topology=md.Topology.from_openmm(simulation.topology)
        python_expression=topology.select_expression(XTC_information)
        req_indices=np.array(eval(python_expression))
        simulation.reporters.append(XTCReporter(outfname, outfreq, atomSubset=req_indices))
    
    ## RUN the simulation
    
    simulation.step(nsteps)
    
    
    ## Output PDB and simulation state
    
    pdb_positions = simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()
    PDBFile.writeFile(simulation.topology, pdb_positions, open(outfile, 'w'))
    simulation.saveState(xmlfile)
    simulation.saveCheckpoint('chkptfile.chk')
    
    ## Get the positions and velocities that could be used to continue the simulation
    
    positions = simulation.context.getState(getPositions=True).getPositions()        #Note PBC condition note enforced -depends on long simulations
    velocities = simulation.context.getState(getVelocities=True).getVelocities()
    
    
    ## Remove the forces added to the simulation by plumed and monte carlo -- better solution will be updated (required for multiple runs in a single script)
    
    if ens=='NPT':
        if use_plumed:
            simulation.context.getSystem().removeForce(simulation.context.getSystem().getNumForces()-1)  #from plumedforce
            simulation.context.getSystem().removeForce(simulation.context.getSystem().getNumForces()-1)  #from Pressure coupl
        else:
            simulation.context.getSystem().removeForce(simulation.context.getSystem().getNumForces()-1)  #from Pressure coupl
            
    return positions,velocities,simulation


############################ Post-hoc saving ############################

def crop_xtc_file(top_file,xtc_file='None',expression='not water',stride=None,skip=0):
    t1=time.perf_counter()
    pdb=md.load(top_file)
    topology=pdb.topology
    python_exp=topology.select_expression(expression)
    req_indices=np.array(eval(python_exp))
    
    top_nosol_file=f"{top_file.split('.')[0]}_nosol.pdb"
    
    md.load(top_file,atom_indices=req_indices).save_pdb(top_nosol_file)
    if xtc_file != 'None':
        xtc_nosol_file=f"{xtc_file.split('.')[0]}_nosol.xtc"
        for i, chunk in enumerate(md.iterload(xtc_file,top=top_file,chunk=1000,atom_indices=req_indices,skip=skip,stride=stride)):
            if i==0:
                trj_nosol=chunk
            else:
                trj_nosol=trj_nosol.join(chunk)
    
        trj_nosol.save_xtc(xtc_nosol_file)
    t2=time.perf_counter()
    print(f'{(t2-t1)/60} mins')

    
