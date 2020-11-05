import scanpy as sc
import pandas as pd
import numpy as np
from pyMAST import MAST

countMatrix = "matrix.mtx"
genesMap = "genes.tsv"
barcodes = "barcodes.tsv"
directory = "hg19"

pathToCountMatrix = f"{directory}/{countMatrix}"
pathToGenesMap = f"{directory}/{genesMap}"
pathToBarcodes = f"{directory}/{barcodes}"


cMatrix = sc.read(pathToCountMatrix, cache=True).T
cMatrix.var_names = pd.read_csv(pathToGenesMap, header=None, sep='\t')[1]
cMatrix.obs_names = pd.read_csv(pathToBarcodes, header=None, sep='\t')[0]

cMatrix.var_names_make_unique()
print(cMatrix)
print(type(cMatrix))
obj = MAST(cMatrix)


