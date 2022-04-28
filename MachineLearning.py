from rdkit import Chem
import deepchem as dc
from transformers import RobertaTokenizerFast
import os

"""
fp_featurizer = dc.feat.CircularFingerprint(size=2048)
pdbDir = "data/pdb-files/"
pdbList = os.listdir(pdbDir)
features = fp_featurizer.featurize([Chem.MolFromPDBFile(pdbDir + pdb) for pdb in pdbList])
"""

tokenizer = RobertaTokenizerFast.from_pretrained("seyonec/PubChem10M_SMILES_BPE_450k")
featurizer = dc.feat.RxnFeaturizer(tokenizer, sep_reagent=True)
feats = featurizer.featurize(['CCS(=O`)(=O)Cl.OCCBr>CCN(CC)CC.CCOCC>CCS(=O)(=O)OCCBr'])
print(feats)
