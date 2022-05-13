# IMPORTS!
from Bio import SeqIO
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit import DataStructs
import logging
# SO MANY I NEED A SECOND COMMENT HALF_WAY DOWN!
import torch
import DataQuery
import numpy as np
import TorchClasses
import pickle
import os.path

# Initiate the proper logging formate
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s: %(message)s")

# If the pickled file versions do not exist, then generate them. Otherwise open them from scratch.
if(not (os.path.isfile("data/tensorRXN.pickle") and os.path.isfile("data/tensorRecords.pickle"))):
    # Initiate the DataFrames
    rhea2uniprot = DataQuery.CreateDF("data/rhea2uniprot_sprot.tsv")
    rheaDirections = DataQuery.CreateDF("data/rhea-directions.tsv")

    # Various variable initializations
    # Will need to go into the for loop in the future
    records = []
    rxnDiffFP = []
    max_seq = 0
    count = 0
    vocab = set()

    # For every record in the given FASTA file,
    for record in SeqIO.parse("data/uniprot-aaseq.fasta", "fasta"):
        # Gather the maximum sequence length
        if(len(str(record.seq)) > max_seq):
            max_seq = len(str(record.seq))

        # update the vocabulary of it 
        vocab.update(str(record.seq))

        # Gather Uniprot and Rhea ID. Without the Rhea ID the entry is useless,
        # and as such thrown aside.
        uniprotID = (str(record.id[3:9]))
        rheaID = DataQuery.QueryDF(uniprotID, "ID", rhea2uniprot)
        if(rheaID.empty):
            continue

        # For every matching Uniprot ID in the record,
        for index, row in rheaID.iterrows():
            # Gather the matching master ID
            rheaRxnID = DataQuery.QueryDF(row["RHEA_ID"], "RHEA_ID_MASTER", rheaDirections)
            # If it's not found, then skipperooni
            if(rheaRxnID.empty):
                continue

            # Otherwise, start checking the direction of the reaction and gather
            # the corresponding LR or RL reaction as necessary - assumed that UN
            # reactions are catalyzed in the LR direction
            elif(row["DIRECTION"] == "UN" or row["DIRECTION"] == "LR"):
                reactionID = rheaRxnID["RHEA_ID_LR"].values[0]
            else:
                reactionID = rheaRxnID["RHEA_ID_RL"]
            
            # Get the reaction from the rxn file into the rdkit Reaction class.
            reaction = AllChem.ReactionFromRxnFile("data/rxn/{}.rxn".format(reactionID))
            
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

            # Append the record and fingerprint to the list.
            records.append(record)
            rxnDiffFP.append(arr)

        # Add one to the count, check progress
        count += 1
        if count % 20 == 0:
            logging.info("{}% finished".format(round(count/9890*100, 2)))

    # Add a pad value to the vocabulary and get a letter-to-number system.
    vocab.add("<pad>")
    to_ix = {char: i for i, char in enumerate(vocab)}

    # Create tensor lists for each.
    tensorRecords = []
    tensorRxnDiff = []

    # Build the tensor record
    for record in records:
        seqTensor = torch.tensor([to_ix[residue] for residue in record.seq]) 
        tensorRecords.append(seqTensor)

    tensorRecords = torch.nn.utils.pad_sequence(
        tensorRecords,
        batch_first=True,
        padding_value=to_ix["<pad>"]
    )
 
    # Build the rxn record
    for rxn in rxnDiffFP:
        rxnTensor = torch.tensor(rxn) 
        tensorRxnDiff.append(rxnTensor)

    # Pickle each for later usage.
    with open("data/tensorRXN.pickle", "wb") as f:
        pickle.dump(rxnTensor, f)

    with open("data/tensorRecords.pickle", "wb") as f:
        pickle.dump(tensorRxnDiff, f)

# If the pickled files already exist, open them!
else:
    with open("data/tensorRXN.pickle", "rb") as f:
        tensorRxnDiff = pickle.load(f)

    with open("data/tensorRecords.pickle", "rb") as f:
        tensorRecords = pickle.load(f)

# Get the number of examples and print as such.
n_examples = len(tensorRecords)
print(n_examples)

# Initialize training and testing datasets
ds_train, ds_test = torch.utils.data.random_split(
    TorchClasses.BiologicalSequenceDataset(tensorRecords, tensorRxnDiff),
    lengths=[n_examples-(n_examples//4), n_examples//4]
)

# Create the corresponding dataloader for each dataset
train_dl = torch.utils.data.DataLoader(ds_train, batch_size=64, shuffle=True)
test_dl = torch.utils.data.DataLoader(ds_test, batch_size=64, shuffle=True)

# Create the neural network and initialize optimizer
nNet = TorchClasses.NeuralNetwork()
criterion = torch.nn.NLLLoss()
optimizer = torch.optim.SGD(nNet.parameters(), lr=0.01, momentum=0.9)

# Perform learning
epochs = 1
for e in range(epochs):
    running_loss = 0
    for reaction, sequence in train_dl:
        optimizer.zero_grad()
        output = nNet(reaction)

        loss = criterion(output, sequence)
        loss.backward()
        optimizer.step()

        running_loss += loss.item()
    else:
        print(f"Training loss: {running_loss/len(train_dl)}")
