import pandas as pd

def CreateDF(tsvFileString):
    return(pd.read_csv(tsvFileString, sep="\t"))

def QueryDF(queryString, queryRow, dfToSearch):
    returnDF = dfToSearch.loc[dfToSearch[queryRow] == queryString]
    return(returnDF)