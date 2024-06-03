import os
import sys
import numpy as np
import pandas as pd
import scipy
import math
import itertools
import pickle

from meeko import MoleculePreparation
from vina import Vina
import subprocess as sp
import mdtraj as md
from rdkit import Chem
from openbabel import openbabel

def mk_lig_pdbqt(lig_file, pdbqt_filepath):
    """ Prepares the ligand pdbqt using Meeko. lig_file is the file path to the pdb. returns the filepath to the new pdbqt. """
    sp.run(['mk_prepare_ligand.py', '-i', lig_file, '-o' ,pdbqt_filepath, '--merge_these_atom_types', '--rigid_macrocycles'])
    return pdbqt_filepath

def protein_pdbqt_ADFR(pdb, pdbqt):
    """ Prepare a protein pdbqt from a pdb using the ADFR suite. """
    sp.run(['ADFR/bin/prepare_receptor' ,'-r', pdb, '-A', 'hydrogens', '-o', pdbqt])
    return pdbqt
  
pdb='/dartfs-hpc/rc/lab/R/RobustelliP/Apara/asn/biorxiv2021-6626290-no-water-glue/lig47.pdb'
trajectory='/dartfs-hpc/rc/lab/R/RobustelliP/Apara/asn/biorxiv2021-6626290-no-water-glue/ligand_47_1.xtc'
dmat = np.load("distance_matrix_full_LIG.npy")
print(dmat.shape)

def combined_residence_events(data):
  idx = np.arange(len(data[:,0]))
  
  bool = np.where(data<.6,1,0)
  any_event_residue = np.any(bool, axis=1)
  
  comp = np.stack([idx,any_event_residue],axis = 1)
  return [i[:,0][1:] if len(i)>1 else i[:,0] for i in filter(lambda x:any(x[:,1]!=0),np.split(comp,np.where(comp[:,1]==0)[0]))]


############################ Start Processing Events ###################

events_all = combined_residence_events(dmat)
print(events_all[:1])

events_all = combined_residence_events(dmat)
mapping_all = np.zeros(len(events_all))
representative_frames_all = np.zeros(len(events_all), dtype=int)
avg = np.zeros([len(events_all), len(dmat[0,:])])
mapping = np.zeros(len(events_all))
representative_frames = np.zeros(len(events_all))
print("Total Events:",len(events_all))


count = 0
for event in events_all:
  mapping_all[count] = len(event)
  representative_frames_all[count] = int(event[int(len(event)/2)])
  count += 1
  
mapping = mapping_all[mapping_all > 4]
representative_frames = representative_frames_all[mapping_all > 4]

print("Total Events At Least 5 Frames:",len(representative_frames))


############################ Getting Docking Score ###################
    
scores = np.zeros(len(mapping))
# scores =  np.load('scores.npy')

for i in range(len(representative_frames)): 
  if scores[i] != 0:
    continue
  frame = representative_frames[i]
  lpath = '/dartfs-hpc/rc/lab/R/RobustelliP/JackM/pdb/l-' + str(frame) + '.pdb'
  lpathmol = '/dartfs-hpc/rc/lab/R/RobustelliP/JackM/mol2/l-' + str(frame) + '.mol2'
  ppath = '/dartfs-hpc/rc/lab/R/RobustelliP/JackM/pdb/p-' + str(frame) + '.pdb'
  
  ligand_file = lpathmol

  protein_file = ppath

  ligand_pdbqt = 'pdbqt_files/l.pdbqt'
  protein_pdbqt = 'pdbqt_files/p.pdbqt'

  mk_lig_pdbqt(ligand_file, ligand_pdbqt)
  protein_pdbqt_ADFR(protein_file, protein_pdbqt)

  try:
    v = Vina(sf_name='vina')
    v.set_receptor('pdbqt_files/p.pdbqt')
    v.set_ligand_from_file('pdbqt_files/l.pdbqt')
    v.compute_vina_maps(center=[0., 0., 0.], box_size=[70, 70, 70])
    scores[i] = np.min(v.score())
  except:
    scores[i] = 999 # flag to remove frame in postprocessing
  
  if i % 50 == 0:
    np.save('scores.npy', scores)
    print(i, mapping[i], scores[i])
    
np.save('scores.npy', scores)