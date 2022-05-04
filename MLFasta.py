from Bio import SeqIO
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit import DataStructs
import logging
import torch
import DataQuery
import numpy as np

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s: %(message)s")


class BiologicalSequenceDataset:
    def __init__(self, sequenceData, rxnData):
        self.seqData = sequenceData
        self.rxnData = rxnData
        self.vocab = set()
        for record in self.seqData:
            self.vocab.update(str(record.seq))
        self.vocab.add("<pad>")
        self.to_ix = {char: i for i, char in enumerate(self.vocab)}


    
    def __len__(self):
        return len(self.seqData)
    
    def __getitem__(self, i):
        seq = torch.tensor([self.to_ix[residue] for residue in self.seqData[i].seq])
        rxn = torch.tensor(self.rxnData[i])
        return seq, rxn

def collate_fn(batch):
    return(torch.nn.utils.rnn.pad_sequence(
        batch,
        batch_first=True,
        padding_value=to_ix["<pad>"]
    ))


rhea2uniprot = DataQuery.CreateDF("data/rhea2uniprot_sprot.tsv")
rheaDirections = DataQuery.CreateDF("data/rhea-directions.tsv")

records = []
rxnDiffFP = []
count = 0

for record in SeqIO.parse("data/uniprot-aaseq.fasta", "fasta"):
    uniprotID = (str(record.id[3:9]))
    rheaID = DataQuery.QueryDF(uniprotID, "ID", rhea2uniprot)

    if(rheaID.empty):
        continue

    for index, row in rheaID.iterrows():
        rheaRxnID = DataQuery.QueryDF(row["RHEA_ID"], "RHEA_ID_MASTER", rheaDirections)
        if(rheaRxnID.empty):
            continue
        elif(row["DIRECTION"] == "UN" or row["DIRECTION"] == "LR"):
            reactionID = rheaRxnID["RHEA_ID_LR"].values[0]
        else:
            reactionID = rheaRxnID["RHEA_ID_RL"]
        
        reaction = AllChem.ReactionFromRxnFile("data/rxn/{}.rxn".format(reactionID))
        
        if(reaction is None):
            continue

        arr = np.zeros((0,), dtype=np.int8)
        rxnFP = rdChemReactions.CreateDifferenceFingerprintForReaction(reaction)
        rxnDiff = DataStructs.ExplicitBitVect(2048)

        for bit in rxnFP:
            rxnDiff.SetBit(bit%rxnDiff.GetNumBits())
        
        DataStructs.ConvertToNumpyArray(rxnDiff, arr)

        if(rxnDiff is None):
            logging.info("BAD. Continuing")
            continue

    count += 1
    if count % 20 == 0:
        logging.info("{}% finished".format(round(count/9890*100, 2)))


n_examples = len(records)
print(n_examples)
ds_train, ds_test = torch.utils.data.random_split(
    BiologicalSequenceDataset(records, rxnDiffFP),
    lengths=[n_examples-(n_examples//4), n_examples//4]
)

train_dl = torch.utils.data.DataLoader(ds_train, batch_size=64, shuffle=True)
test_dl = torch.utils.data.DataLoader(ds_test, batch_size=64, shuffle=True)

