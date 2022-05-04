from Bio import SeqIO
from rdkit import Chem
from rdkit.Chem import AllChem
import deepchem as dc
import numpy as np
import os
import logging
from transformers import RobertaTokenizerFast

import DataQuery

print()
print()
rhea2uniprot = DataQuery.CreateDF("data/rhea2uniprot_sprot.tsv")
rheaDirections = DataQuery.CreateDF("data/rhea-directions.tsv")
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s: %(message)s")

logging.info('Beginning script')
count = 0
codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
seqList = []
rxnSmartsList = []
uniprotIDs = []
maxseq = 0

logging.info("Opening seqarr and rxnarr")
if not(os.path.isfile("data/seqarr.npy") and os.path.isfile("data/rxnarr.npy")):
    logging.info("One is not found, featurizing from scratch")
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

            reactionSmarts = Chem.rdChemReactions.ReactionToSmarts(reaction)

            if(reactionSmarts is None or str(reactionSmarts).strip() == ""):
                print("BAD. Continuing")
                continue

            seqList.append(str(record.seq))
            rxnSmartsList.append(reactionSmarts)
            uniprotIDs.append(uniprotID)

        if(maxseq < len(record.seq)):
            maxseq = len(record.seq)
        
        count += 1
        if(count % 10 == 0):
            print("{}% complete".format(round(count/9890*100, 2)), end='\r')
        
    logging.info("Finished loading data, now featurizing or opening pre-featurized...")
    logging.info("The maximum sequence length is {}".format(maxseq))
    if(not os.path.isfile("data/seqarr.npy")):
        logging.info("Seq List Length: {}".format(len(seqList)))
        aafeat = dc.feat.OneHotFeaturizer(codes, max_length=maxseq)
        logging.info("Featurizing sequences...")
        output = aafeat.featurize(seqList[0:9000])
        logging.info("Saving the sequences to data/seqarr.npy")
        with open("data/seqarr.npy", "wb") as f:
            np.save(f, output)
    else:
        with open("data/seqarr.npy", "rb") as f:
            output = np.load(f)

    print("Finished featurizing sequences, now for reactions...")
    if(not os.path.isfile("data/rxnarr.npy")):
        tokenizer = RobertaTokenizerFast.from_pretrained("seyonec/PubChem10M_SMILES_BPE_450k")
        rxnfeat = dc.feat.RxnFeaturizer(tokenizer, sep_reagent=True)
        features = rxnfeat.featurize(rxnSmartsList[0:9000])
        with open("data/rxnarr.npy", "wb") as f:
            np.save(f, features)
    else:
        with open("data/rxnarr.npy", "rb") as f:
            features = np.load(f)
else:
    with open("data/rxnarr.npy", "rb") as f:
        features = np.load(f, allow_pickle=True)
    with open("data/seqarr.npy", "rb") as f:
        output = np.load(f, allow_pickle=True)
        
logging.info("Creating Deepchem Dataset from created data")
dataset = dc.data.DiskDataset.from_numpy(X=features, y=output, data_dir="data/dataset")