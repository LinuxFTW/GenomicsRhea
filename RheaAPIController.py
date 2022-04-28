import requests
import pandas as pd
import xml.etree.ElementTree as ET
from rdflib import Graph

rheaURL = "https://www.rhea-db.org/rhea?"
def getRXN(rheaID):
    params = {
        "query":'',
        "columns":"rhea-id",
        "format":'tsv',
        "limit":10
    }

def readRDF(rdfFile):
    tree = ET.ElementTree(file=rdfFile)
    print("Finished processing in XML, now it's RDF")
    tree2 = Graph()
    tree2.parse("data/rhea.rdf")
    for term in tree2:
        print(term)

readRDF("data/rhea.rdf")

