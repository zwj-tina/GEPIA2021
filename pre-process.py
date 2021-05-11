import pandas as pd
import numpy as np

# tcga
data = pd.read_csv("./tcga_RSEM_gene_tpm", sep = "\t")# gtex 
gene_map = pd.read_csv("./probeMap%2Fgencode.v23.annotation.gene.probemap", sep = "\t")
ids = gene_map["id"].tolist()
gene = gene_map["gene"].tolist()
geneSymbol = []
for each in data["sample"].tolist():
    po = ids.index(each)
    geneSymbol.append(gene[po])
data["geneSymbol"] = geneSymbol
data = data.drop_duplicates(['geneSymbol'])
indexes = data["geneSymbol"].tolist()
data = data.iloc[:,1:-1]
data.index = indexes
data_tpm = 2 ** data  - 0.001
protein =  pd.read_csv("./protein.txt", sep = "\t")
gene = protein["Symbol"].tolist()
gene_need = [i for i in gene if i in data_tpm.index.tolist()]
data_tpm = data_tpm.loc[gene_need,:]
col_sum = data_tpm.sum(axis = 0)
ax = 1e6 / col_sum
data_tpm = data_tpm * ax
data_tpm.to_csv("./tcgaTpm.csv",index_label = "gene")


#gtex
data = pd.read_csv("./gtex_RSEM_gene_tpm", sep = "\t")# gtex 
ids = gene_map["id"].tolist()
gene = gene_map["gene"].tolist()
geneSymbol = []
for each in data["sample"].tolist():
    po = ids.index(each)
    geneSymbol.append(gene[po])
data["geneSymbol"] = geneSymbol
data = data.drop_duplicates(['geneSymbol'])
indexes = data["geneSymbol"].tolist()
data = data.iloc[:,1:-1]
data.index = indexes
data_tpm = 2 ** data  - 0.001
protein =  pd.read_csv("./protein.txt", sep = "\t")
gene = protein["Symbol"].tolist()
gene_need = [i for i in gene if i in data_tpm.index.tolist()]
data_tpm = data_tpm.loc[gene_need,:]
col_sum = data_tpm.sum(axis = 0)
ax = 1e6 / col_sum
data_tpm = data_tpm * ax
data_tpm.to_csv("./gtexTpm.csv",index_label = "gene")
