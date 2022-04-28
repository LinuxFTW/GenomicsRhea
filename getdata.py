# The purpose of this file is data extraction - getting all of the PDB and other
# data required to do the machine learning.
import pandas as pd
import re
import urllib.request
import urllib.parse
import os
import wget
import gzip
import shutil
from ReverseQuery import *

print("Initializing variables and loading databases")
# Variable initialization - for pdb, rheadb, and the url including www. because apparently that's important.
# Location of PDBs
pdbDir = "data/pdb-files/"

# List of currently available PDBs in pdbDir
availablePDBs = os.listdir(pdbDir)

# Rhea database loaded into Pandas DF format and Rhea2Uniprot database loaded into
# Pandas DF.
rheadb = pd.read_csv("data/rhea-ec-iubmb.tsv", sep="\t")
rhea2uniprot = pd.read_csv("data/rhea2uniprot_sprot.tsv", sep="\t")

# Regex string
ecMatch = "^1."

# API links
uniprotAPI = 'https://www.uniprot.org/uploadlists/'
pdbAPI = 'https://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/'

# List of .rxn files
rxnList = os.listdir("data/rd")

# Checking to see if uniprot-pdb has already been made, and if not, gets all of the data to create it
if("uniprot-pdb.tsv" not in os.listdir("data/")):
    print("Uniprot-PDB Database not found, generating...")
    uniprotIDstring = ''
    # Iterate through each row in rheadb, and if they are a hydrolysis reaction
    # get its associated uniprotID and add it to the query list.
    for index, row in rheadb.iterrows():
        print(row["REACTION_ID"], str(row["REACTION_ID"]) + ".rd" in rxnList)
        if(re.search(ecMatch, row["EC"])):
            uniprotID = QueryDF(row["REACTION_ID"], "RHEA_ID", rhea2uniprot)
            print(uniprotID)

            if(uniprotID.empty):
                continue
            
            for index2, enzyme in uniprotID.iterrows():
                if(enzyme["ID"] not in availablePDBs and str(row["REACTION_ID"]) + ".rxn" in rxnList):
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
    req = urllib.request.Request(uniprotAPI, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    
    # Write the data
    with open("data/uniprot-pdb.tsv", "w") as uniprotPDB:
        uniprotPDB.write(response.decode('utf-8'))
else:
    print("Found Uniprot-PDB Databse...")

print("Downloading PDBS not found")
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
    else:
        print("{} found in {}".format(pdbEntry, pdbDir))
