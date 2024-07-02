import numpy as np
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import Chem
from rxnmapper import RXNMapper

reactants_patt_list=['[F,Cl,Br,I]-[c;H0;D3;+0:1](:[c,n:2]):[c,n:3].O-B(-O)-[c;H0;D3;+0:4](:[c,n:5]):[c,n:6]']
product_patt_list=['[c,n:2]:[c;H0;D3;+0:1](:[c,n:3])-[c;H0;D3;+0:4](:[c,n:5]):[c,n:6]']
reactants_protect_patt_list=['[F,Cl,Br,I]-[c;H0;D3;+0:1](:[c,n:2]):[c,n:3]-C1OCCO1.O-B(-O)-[c;H0;D3;+0:4](:[c,n:5]):[c,n:6]']
product_protect_patt_list=['O1CCOC1-[c,n:2]:[c;H0;D3;+0:1](:[c,n:3])-[c;H0;D3;+0:4](:[c,n:5]):[c,n:6]']
reactants_without_protect_patt_list=['[F,Cl,Br,I]-[c;H0;D3;+0:1](:[c,n:2]):[c,n:3]-C=O.O-B(-O)-[c;H0;D3;+0:4](:[c,n:5]):[c,n:6]']
product_without_protect_patt_list=['O=C-[c,n:2]:[c;H0;D3;+0:1](:[c,n:3])-[c;H0;D3;+0:4](:[c,n:5]):[c,n:6]']
reactants_NH2_without_protect_patt_list=['[F,Cl,Br,I]-[c;H0;D3;+0:1](:[c,n:2]):[c,n:3]-N.O-B(-O)-[c;H0;D3;+0:4](:[c,n:5]):[c,n:6]']
product_NH2_without_protect_patt_list=['N-[c,n:2]:[c;H0;D3;+0:1](:[c,n:3])-[c;H0;D3;+0:4](:[c,n:5]):[c,n:6]']
reactants_OH_without_protect_patt_list=['[F,Cl,Br,I]-[c;H0;D3;+0:1](:[c,n:2]):[c,n:3]-O.O-B(-O)-[c;H0;D3;+0:4](:[c,n:5]):[c,n:6]']
product_OH_without_protect_patt_list=['O-[c,n:2]:[c;H0;D3;+0:1](:[c,n:3])-[c;H0;D3;+0:4](:[c,n:5]):[c,n:6]']


def find_substructure(m,patt_list):
    for patt in patt_list:
        patt = Chem.MolFromSmarts(patt)
        flag = m.HasSubstructMatch(patt)
        if flag:
            break
    return flag


def find_suzuki_template(rxn):
    #print(rxn)
    if rxn.count('>>')==1:    
        if rxn.count('.')==1:
            reactants, products = rxn.split('>>')
            if "." in reactants:
                reactants_mol=Chem.MolFromSmiles(reactants)
                if reactants_mol is None:
                    reactants_find_result=False
                else:
                    reactants_find_result=find_substructure(reactants_mol, reactants_patt_list)
                    
                product_mol=Chem.MolFromSmiles(products)
                if product_mol is None:
                    product_find_result=False
                else:
                    product_find_result=find_substructure(product_mol,product_patt_list)
            else:
                reactants_find_result=False
                product_find_result=False
        else:
            reactants_find_result=False
            product_find_result=False
    else:
        reactants_find_result=False
        product_find_result=False
    final_result=reactants_find_result&product_find_result
    return final_result  

def find_suzuki_protect_template(rxn):
    #print(rxn)
    if rxn.count('>>')==1:    
        if rxn.count('.')==1:
            reactants, products = rxn.split('>>')
            if "." in reactants:
                reactants_mol=Chem.MolFromSmiles(reactants)
                if reactants_mol is None:
                    reactants_find_result=False
                else:
                    reactants_find_result=find_substructure(reactants_mol, reactants_protect_patt_list)
                    
                product_mol=Chem.MolFromSmiles(products)
                if product_mol is None:
                    product_find_result=False
                else:
                    product_find_result=find_substructure(product_mol,product_protect_patt_list)
            else:
                reactants_find_result=False
                product_find_result=False
        else:
            reactants_find_result=False
            product_find_result=False
    else:
        reactants_find_result=False
        product_find_result=False
    final_result=reactants_find_result&product_find_result
    return final_result

def find_suzuki_without_protect_template(rxn):
    #print(rxn)
    if rxn.count('>>')==1:    
        if rxn.count('.')==1:
            reactants, products = rxn.split('>>')
            if "." in reactants:
                reactants_mol=Chem.MolFromSmiles(reactants)
                if reactants_mol is None:
                    reactants_find_result=False
                else:
                    reactants_find_result=find_substructure(reactants_mol, reactants_without_protect_patt_list)
                    
                product_mol=Chem.MolFromSmiles(products)
                if product_mol is None:
                    product_find_result=False
                else:
                    product_find_result=find_substructure(product_mol,product_without_protect_patt_list)
            else:
                reactants_find_result=False
                product_find_result=False
        else:
            reactants_find_result=False
            product_find_result=False
    else:
        reactants_find_result=False
        product_find_result=False
    final_result=reactants_find_result&product_find_result
    return final_result


def find_suzuki_NH2_without_protect_template(rxn):
    #print(rxn)
    if rxn.count('>>')==1:    
        if rxn.count('.')==1:
            reactants, products = rxn.split('>>')
            if "." in reactants:
                reactants_mol=Chem.MolFromSmiles(reactants)
                if reactants_mol is None:
                    reactants_find_result=False
                else:
                    reactants_find_result=find_substructure(reactants_mol, reactants_NH2_without_protect_patt_list)
                    
                product_mol=Chem.MolFromSmiles(products)
                if product_mol is None:
                    product_find_result=False
                else:
                    product_find_result=find_substructure(product_mol,product_NH2_without_protect_patt_list)
            else:
                reactants_find_result=False
                product_find_result=False
        else:
            reactants_find_result=False
            product_find_result=False
    else:
        reactants_find_result=False
        product_find_result=False
    final_result=reactants_find_result&product_find_result
    return final_result

def find_suzuki_OH_without_protect_template(rxn):
    #print(rxn)
    if rxn.count('>>')==1:    
        if rxn.count('.')==1:
            reactants, products = rxn.split('>>')
            if "." in reactants:
                reactants_mol=Chem.MolFromSmiles(reactants)
                if reactants_mol is None:
                    reactants_find_result=False
                else:
                    reactants_find_result=find_substructure(reactants_mol, reactants_OH_without_protect_patt_list)
                    
                product_mol=Chem.MolFromSmiles(products)
                if product_mol is None:
                    product_find_result=False
                else:
                    product_find_result=find_substructure(product_mol,product_OH_without_protect_patt_list)
            else:
                reactants_find_result=False
                product_find_result=False
        else:
            reactants_find_result=False
            product_find_result=False
    else:
        reactants_find_result=False
        product_find_result=False
    final_result=reactants_find_result&product_find_result
    return final_result