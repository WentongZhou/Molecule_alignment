import copy
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
from rdkit.Chem import rdBase
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolDescriptors
import numpy as np
p = AllChem.ETKDGv3()
p.verbose = True
import argparse

def align_mols(ref_mol,mols,method='crippen'):
    ref_mol = Chem.SDMolSupplier(ref_mol)[0]
    ref_mol = Chem.AddHs(ref_mol)
    AllChem.EmbedMultipleConfs(ref_mol, 100, p)
    mols = [m for m in Chem.SDMolSupplier(mols) if m != None]
    for mol in mols:
        mol.RemoveAllConformers()
    hmols_1 = [Chem.AddHs(m) for m in mols]
    hmols_2 = copy.deepcopy(hmols_1)
    for mol in hmols_1:
        AllChem.EmbedMultipleConfs(mol, 100, p)
    for mol in hmols_2:
        AllChem.EmbedMultipleConfs(mol, 100, p)

    Chem.MolToXYZFile(ref_mol, f'ref.xyz')
    if method == 'crippen':
        crippen_contribs = [rdMolDescriptors._CalcCrippenContribs(mol) for mol in hmols_1]
        crippen_ref_contrib = rdMolDescriptors._CalcCrippenContribs(ref_mol)
        crippen_prob_contribs = crippen_contribs
        prob_mols_1 = hmols_1
        for idx, mol in enumerate(prob_mols_1):
            tempscore = []
            for cid in range(100):
                crippenO3A = rdMolAlign.GetCrippenO3A(mol, ref_mol, crippen_prob_contribs[idx], crippen_ref_contrib, cid, 0)
                crippenO3A.Align()
                tempscore.append(crippenO3A.Score())
            best = np.argmax(tempscore)
            Chem.MolToXYZFile(mol, confId=int(best), filename=f'{idx}.xyz')
    elif method == 'mmff':
        mmff_params = [AllChem.MMFFGetMoleculeProperties(mol) for mol in hmols_2]
        mmff_ref_param = AllChem.MMFFGetMoleculeProperties(ref_mol)
        mmff_prob_params = mmff_params
        prob_mols_2 = hmols_2
        for idx, mol in enumerate(prob_mols_2):
            tempscore = []
            for cid in range(100):
                pyO3A = rdMolAlign.GetO3A(mol, ref_mol, mmff_prob_params[idx], mmff_ref_param, cid, 0)
                pyO3A.Align()
                tempscore.append(pyO3A.Score())
            best = np.argmax(tempscore)
            Chem.MolToXYZFile(mol, confId=int(best), filename=f'{idx}.xyz')

if __name__ == '__main__':
    # generate arguments
    parser = argparse.ArgumentParser(description='Align molecules')
    parser.add_argument('-r', '--ref', type=str, help='reference molecule')
    parser.add_argument('-m', '--mols', type=str, help='molecules to be aligned')
    parser.add_argument('-t', '--type', type=str, help='alignment method')
    args = parser.parse_args()
    # run the function with arguments
    align_mols(args.ref,args.mols,args.type)
