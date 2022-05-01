from rdkit import Chem
import deepchem as dc
from transformers import RobertaTokenizerFast
import os
import time

fp_featurizer = dc.feat.ContactCircularFingerprint(size=2048)
pdbDir = "data/pdb-files/"
pdbList = os.listdir(pdbDir)
startTime = time.time()
MolsFromPDB = [Chem.MolFromPDBFile(pdbDir+pdb) for pdb in pdbList]
features = fp_featurizer.featurize(pdbList)
stopTime = time.time()

print(features)
print("It took {} seconds to execute".format(stopTime-startTime))

"""
tokenizer = RobertaTokenizerFast.from_pretrained("seyonec/PubChem10M_SMILES_BPE_450k")
featurizer = dc.feat.RxnFeaturizer(tokenizer, sep_reagent=True)
feats = featurizer.featurize(['CCS(=O`)(=O)Cl.OCCBr>CCN(CC)CC.CCOCC>CCS(=O)(=O)OCCBr'])
print(feats)
"""


