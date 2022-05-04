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
import ReactionTest
from DataQuery import *

print("Initializing variables and loading databases")

# Location of PDBs
pdbDir = "data/pdb-files/"

# List of currently available PDBs in pdbDir
availablePDBs = os.listdir(pdbDir)

# Rhea database loaded into Pandas DF format and Rhea2Uniprot database loaded into
# Pandas DF.
rhea2ec = CreateDF("data/rhea-ec-iubmb.tsv")
rhea2uniprot = CreateDF("data/rhea2uniprot_sprot.tsv")
rhea2rxndir = CreateDF("data/rhea-directions.tsv")

# Regex string
# HYDROLASES ARE 3. not 1. AAAAAAAAAA
ecMatch = "^3."

# API links
uniprotAPI = 'https://www.uniprot.org/uploadlists/'
pdbAPI = 'https://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/'

# List of .rxn files
rxnList = os.listdir("data/rxn")

# Checking to see if uniprot-pdb has already been made, and if not, gets all of the data to create it
if("uniprot-pdb.tsv" not in os.listdir("data/")):
    print("Uniprot-PDB Database not found, generating...")
    uniprotIDstring = ''
    unavailableReactions = ReactionTest.CheckReactions()
    # Iterate through each row in rhea2ec, and if they are a hydrolysis reaction
    # get its associated uniprotID and add it to the query list.
    for index, row in rhea2ec.iterrows():
        if(re.search(ecMatch, row["EC"])):
            uniprotID = QueryDF(row["REACTION_ID"], "RHEA_ID", rhea2uniprot)
            L2RDirection = QueryDF(row["REACTION_ID"], "RHEA_ID_MASTER", rhea2rxndir)

            if(uniprotID.empty or L2RDirection.empty):
                continue

            for index2, enzyme in uniprotID.iterrows():
                uniprotIDstring += enzyme["ID"] + ' '
        if(index % 10):
            print("This process is {}% done".format(round(index/len(rhea2ec) * 100, 2)), end="\r")
    
    print("Querying Uniprot for available data...")
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

# Download unfound PDBs
print("Downloading PDBS not found")
uniprotPDB = CreateDF("data/uniprot-pdb.tsv")

for index, row in uniprotPDB.iterrows():
    # Generate the corresponding PDB entry
    pdbEntry = "pdb" + row["To"].lower() + ".ent.gz"
    # If said entry is not in the list, generate the link and download.
    if(pdbEntry[:-3] not in availablePDBs):
        pdbDL = pdbAPI + pdbEntry[4:6] + "/" + pdbEntry
        try:
            download = wget.download(pdbDL, out=pdbDir)
        except urllib.error.HTTPError:
            continue
        except urllib.error.URLError:
            continue
        
        # Decompress downloaded .ent.gz file and write it out. Remove
        # the .ent file at the end of the download.
        with gzip.open(pdbDir + "/" + pdbEntry) as pdbCompressed:
            with open(pdbDir + pdbEntry[:-3], 'wb') as pdbDecompressed:
                shutil.copyfileobj(pdbCompressed, pdbDecompressed)
        os.remove(pdbDir + pdbEntry)
    else:
        print("{} found in {}".format(pdbEntry, pdbDir))