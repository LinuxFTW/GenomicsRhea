# The purpose of this file is data extraction - getting all of the PDB and other
# data required to do the machine learning.
from rdkit import Chem
import deepchem as dc
import tensorflow as tf
import pandas as pd
import re
import urllib.request
import urllib.parse
import os
import wget
import gzip
import shutil

# Variable initialization - for pdb, rheadb, and the url including www. because apparently that's important.
pdbDir = "data/pdb-files/"
availablePDBs = os.listdir(pdbDir)
rheadb = pd.read_csv("data/rhea-ec-iubmb.tsv", sep="\t")
rhea2uniprot = pd.read_csv("data/rhea2uniprot_sprot.tsv", sep="\t")
ecMatch = "^1."
uniprotAPI = 'https://www.uniprot.org/uploadlists/'
pdbAPI = 'https://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/'

# Checking to see if uniprot-pdb has already been made, and if not, gets all of the data to create it
print(os.listdir("data/"))
if("uniprot-pdb.tsv" not in os.listdir("data/")):
    print("Uniprot-PDB Database not found, generating...")
    uniprotIDstring = ''
    # Iterate through each row in rheadb, and if they are a hydrolysis reaction
    # get its associated uniprotID and add it to the query list.
    for index, row in rheadb.iterrows():
        if(re.search(ecMatch, row["EC"])):
            uniprotID = rhea2uniprot.loc[rhea2uniprot["RHEA_ID"] == row["REACTION_ID"]]

            if(uniprotID.empty):
                continue
            
            for index2, enzyme in uniprotID.iterrows():
                if(enzyme["ID"] not in availablePDBs):
                    uniprotIDstring += enzyme["ID"] + ' '
    
    # Initialize params
    params = {
    'from': 'ACC+ID',
    'to': 'PDB_ID',
    'format': 'tab',
    'query': uniprotIDstring
    }

    # Get the data
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    
    # Write the data
    with open("data/uniprot-pdb.tsv", "w") as uniprotPDB:
        uniprotPDB.write(response.decode('utf-8'))

uniprotPDB = pd.read_csv("data/uniprot-pdb.tsv", delimiter="\t")
for index, row in uniprotPDB.iterrows():
    pdbEntry = "pdb" + row["To"].lower() + ".ent.gz"
    if(pdbEntry[:-3] not in availablePDBs):
        pdbDL = pdbAPI + pdbEntry[4:6] + "/" + pdbEntry
        try:
            download = wget.download(pdbDL, out=pdbDir)
        except urllib.error.HTTPError:
            continue
        except urllib.error.URLError:
            continue

        with gzip.open(pdbDir + "/" + pdbEntry) as pdbCompressed:
            with open(pdbDir + pdbEntry[:-3], 'wb') as pdbDecompressed:
                shutil.copyfileobj(pdbCompressed, pdbDecompressed)
        os.remove(pdbDir + pdbEntry)
