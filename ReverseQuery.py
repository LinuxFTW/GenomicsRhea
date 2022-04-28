import pandas as pd

def QueryDF(queryString, queryRow, dfToSearch):
    returnDF = dfToSearch.loc[dfToSearch[queryRow] == queryString]
    return(returnDF)