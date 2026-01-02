import pandas as pd

def namer(sampleTable, columnNames):
    nameList = []
    # For columns which you want to list ONLY if there is a single one
    for column in columnNames:
        columnSet = list(set(sampleTable[column]))
        if len(columnSet) == 1:
            nameList.append(columnSet[0])
    return("_".join(nameList))

def projectnamer(sampleTable, columnNames):
    nameList = []
    for column in columnNames:
        columnSet = sorted(set(sampleTable[column]))
        if len(columnSet) == 1:
            nameList.append(str(columnSet[0]))
        else:
            unique_values = '_'.join(str(value) for value in columnSet if pd.notna(value))
            if unique_values:
                nameList.append(f"{column}-{unique_values}")
    return "_".join(nameList) if nameList else "AllSamples"
