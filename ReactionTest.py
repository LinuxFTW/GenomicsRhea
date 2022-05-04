from rdkit.Chem import AllChem
import os

def CheckReactions():
    rxnFolder = os.listdir("data/rxn/")
    failedReactions = []

    for rxnFile in rxnFolder:
        try:
            reaction = AllChem.ReactionFromRxnFile("data/rxn/" + rxnFile)
        except (RuntimeError, TypeError, ValueError):
            failedReactions.append(rxnFile[:-4])
    return(failedReactions)
