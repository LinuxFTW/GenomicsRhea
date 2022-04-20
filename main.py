# The purpose of this file is data extraction - getting all of the PDB and other
# data required to do the machine learning.
from rdkit import Chem
import deepchem as dc
import tensorflow as tf
import pandas as pd
import re
import requests

dataDir = "data/"
rheadb = pd.read_csv(dataDir+"rhea-ec-iubmb.tsv", sep="\t")
ecMatch = "^1."
rheaAPI = "https://www.rhea-db.org/rhea?"
rheaAPI_params = {
    "query":'uniprot:*',
    "columns":"rhea-id,equation,uniprot",
    "format":"tsv",
    "limit":10
}

response = requests.get(rheaAPI,params=rheaAPI_params)
print(dir(response))
for index, row in rheadb.iterrows():
    continue