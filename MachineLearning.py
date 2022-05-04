from rdkit import Chem
import deepchem as dc
from transformers import RobertaTokenizer
import os
import time
import numpy as np

from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.PDBParser import PDBParser

fp_featurizer = dc.feat.ConvMolFeaturizer(per_atom_fragmentation=False)
pdbDir = "data/pdb-files/"
pdbList = os.listdir(pdbDir)
startTime = time.time()

MolsFromPDB = []
count = 0
total = len(pdbList)

for pdb in pdbList[0:100]:
    print("a")
    structure = PDBParser.get_structure(pdb[3:7], pdbDir+pdb)
    abc = PPBuilder().build_peptides(structure)
    print(abc)

"""
for pdb in pdbList[0:100]:
    rdkitFromPDB = Chem.MolFromPDBFile(pdbDir+pdb)
    if(rdkitFromPDB is None):
        print("Cannot add {} as it did not convert properly".format(pdbDir+pdb))
        continue
    MolsFromPDB.append(rdkitFromPDB)
    count += 1
    if count % 10:
        print("{}% completed".format(count/total*100), end="\r")

print("Done converting proteins to molecules, featurizing")

features = fp_featurizer.featurize([Chem.MolToSmiles(pdb) for pdb in MolsFromPDB])
stopTime = time.time()

print(features)
print("It took {} seconds to execute".format(stopTime-startTime))

tokenizer = RobertaTokenizerFast.from_pretrained("seyonec/PubChem10M_SMILES_BPE_450k")
featurizer = dc.feat.RxnFeaturizer(tokenizer, sep_reagent=True)
feats = featurizer.featurize(['CCS(=O`)(=O)Cl.OCCBr>CCN(CC)CC.CCOCC>CCS(=O)(=O)OCCBr'])
print(feats)

"""
