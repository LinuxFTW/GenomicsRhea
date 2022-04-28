from rdkit.Chem import AllChem
import os

rxnFolder = os.listdir("data/rxn/")

for rxnFile in rxnFolder:
    try:
        reaction = AllChem.ReactionFromRxnFile("data/rxn/" + rxnFile)
        print(AllChem.ReactionToSmarts(reaction))
    except (RuntimeError, TypeError, ValueError):
        print(rxnFile + " failed. Continuing")