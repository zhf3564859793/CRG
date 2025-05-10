import numpy as np
import pandas as pd
import warnings
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdChemReactions import RemoveMappingNumbersFromReactions

def get_rxn_maccs(rxn):
    rxn = AllChem.ReactionFromSmarts(rxn,useSmiles=True)
    rxn1=RemoveMappingNumbersFromReactions(rxn)
    rxn=AllChem.ReactionToSmiles(rxn)
    rt,pt = rxn.split('>>')[0], rxn.split('>>')[1]
    
    mol_rt=Chem.MolFromSmiles(rt)
    mol_pt=Chem.MolFromSmiles(pt)
    rt_maccs_fps=MACCSkeys.GenMACCSKeys(mol_rt)
    pt_maccs_fps=MACCSkeys.GenMACCSKeys(mol_pt)
    
    return [rt_maccs_fps,pt_maccs_fps]

def cal_rxn_dice_similarity(r1,r2):
    sim_r=DataStructs.FingerprintSimilarity(r1[0],r2[0],metric=DataStructs.DiceSimilarity)
    sim_p=DataStructs.FingerprintSimilarity(r1[1],r2[1],metric=DataStructs.DiceSimilarity)
    sim=(sim_r+sim_p)/2
    return sim

def get_max_similarity(rxn,rxn_list):
    sim_list=[]
    for j in rxn_list:
        sim=cal_rxn_dice_similarity(rxn,j)
        sim_list.append(sim)
    max_sim=max(sim_list)
    return max_sim

def get_max_similarity_1(rxn,rxn_list):
    sim_list=[]
    for j in rxn_list:
        sim=cal_rxn_dice_similarity(rxn,j)
        if sim!=1:
            sim_list.append(sim)
    max_sim=max(sim_list)
    return max_sim

def cal_rxn_tanimoto_similarity(r1,r2):
    sim_r=DataStructs.FingerprintSimilarity(r1[0],r2[0],metric=DataStructs.TanimotoSimilarity)
    sim_p=DataStructs.FingerprintSimilarity(r1[1],r2[1],metric=DataStructs.TanimotoSimilarity)
    sim=(sim_r+sim_p)/2
    return sim

def get_max_tanimoto_similarity(rxn,rxn_list):
    sim_list=[]
    for j in rxn_list:
        sim=cal_rxn_tanimoto_similarity(rxn,j)
        sim_list.append(sim)
    max_sim=max(sim_list)
    return max_sim
    

def get_rxn_ecfp(rxn):
    rxn = AllChem.ReactionFromSmarts(rxn,useSmiles=True)
    rxn1=RemoveMappingNumbersFromReactions(rxn)
    rxn=AllChem.ReactionToSmiles(rxn)
    rt,pt = rxn.split('>>')[0], rxn.split('>>')[1]
    
    mol_rt=Chem.MolFromSmiles(rt)
    mol_pt=Chem.MolFromSmiles(pt)
    rt_ecfp_fps=AllChem.GetMorganFingerprintAsBitVect(mol_rt,2,nBits=1024)
    pt_ecfp_fps=AllChem.GetMorganFingerprintAsBitVect(mol_pt,2,nBits=1024)
    
    return [rt_ecfp_fps,pt_ecfp_fps]

def cal_rxn_ecfp_dice_similarity(r1,r2):
    sim_r=DataStructs.DiceSimilarity(r1[0],r2[0])
    sim_p=DataStructs.DiceSimilarity(r1[1],r2[1])
    sim=(sim_r+sim_p)/2
    return sim


def get_max_ecfp_dice_similarity(rxn,rxn_list):
    sim_list=[]
    for j in rxn_list:
        sim=cal_rxn_ecfp_dice_similarity(rxn,j)
        
        sim_list.append(sim)
    max_sim=max(sim_list)
    return max_sim


def cal_rxn_ecfp_tanimoto_similarity(r1,r2):
    sim_r=DataStructs.FingerprintSimilarity(r1[0],r2[0],metric=DataStructs.TanimotoSimilarity)
    sim_p=DataStructs.FingerprintSimilarity(r1[1],r2[1],metric=DataStructs.TanimotoSimilarity)
    sim=(sim_r+sim_p)/2
    return sim

def get_max_ecfp_tanimoto_similarity(rxn,rxn_list):
    sim_list=[]
    for j in rxn_list:
        sim=cal_rxn_ecfp_tanimoto_similarity(rxn,j)    
        sim_list.append(sim)
    max_sim=max(sim_list)
    return max_sim

