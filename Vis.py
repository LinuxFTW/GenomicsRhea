import torch
import TorchClasses
import os
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit import DataStructs
import numpy as np

input_names = ["Reaction Difference"]
output_names = ["Featurized Sequence"]

model = TorchClasses.NeuralNetwork()
model.load_state_dict(torch.load("data/NeuralNet99.pickle"))
screaming = []

for reactionID in os.listdir("data/rxn"):
    print("a")
   # Get the reaction from the rxn file into the rdkit Reaction class.
    try:
       reaction = AllChem.ReactionFromRxnFile("data/rxn/{}".format(reactionID))
    except ValueError:
       continue
    print("b")

   # Ignore it if RDKit does not like the data presentation
    if(reaction is None):
       continue

   # Initialize variables and greate the diference fingerprint for the reaction.
    arr = np.zeros((0,), dtype=np.int8)
    rxnFP = rdChemReactions.CreateDifferenceFingerprintForReaction(reaction)
    rxnDiff = DataStructs.ExplicitBitVect(2048)

   # For every bit in the fingerprint, translate it to the Explicit Bit Vector format.
    for bit in rxnFP:
        rxnDiff.SetBit(bit%rxnDiff.GetNumBits())

    # Convert that Explicit Bit Vector to the Numpy array
    DataStructs.ConvertToNumpyArray(rxnDiff, arr)

    # Ensure rxnDiff is a real variable
    if(rxnDiff is None):
        logging.info("BAD. Continuing")
        continue

    screaming.append(torch.tensor(rxnDiff).float())
    break

unuseful = model(screaming[0])

torch.onnx.export(model, screaming[0], "data/model.onnx", input_names=input_names, output_names=output_names)