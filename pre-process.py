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
data_tpm.to_csv("./tcgaTpmGenes.csv",index_label = "gene")
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
data_tpm.to_csv("./gtexTpmGenes.csv",index_label = "gene")
