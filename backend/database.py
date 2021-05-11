import os
import numpy as np
import pandas as pd
import API

#tcga
tcga_matrix = pd.read_csv("./tcgaTpm.csv")
genes = tcga_matrix["gene"].tolist()
tcga_matrix_t = tcga_matrix.iloc[:,1:].T
tcga_matrix_t.columns = genes
data = pd.read_csv("./TCGA_phenotype_denseDataOnlyDownload.tsv",sep="\t")
propotion = pd.read_csv("./tcga_EPIC_propotion.csv")
cellID = propotion["Unnamed: 0"].tolist()
quan_propotion = pd.read_csv("./tcga_QS_propotion.csv")
cell_type = quan_propotion["cell_type"].tolist()
quan_propotion = quan_propotion.iloc[:,2:].T
quan_propotion.columns = cell_type
LM_propotion = pd.read_csv("./tcga_LM_propotion.csv")
cols = LM_propotion["Unnamed: 0"].tolist()
cell_type = LM_propotion.columns.tolist()[1:]
LM_propotion = LM_propotion.iloc[:,1:]
LM_propotion.columns = cell_type
LM_propotion.index = cols
LM_n = LM_propotion.iloc[:,:22]
sample = data["sample"].tolist()
sample_correct = []
for each in sample:
    sample_correct.append(each.replace("-","."))
data.index = sample_correct
union = [i for i in sample_correct if i in cellID]
data_redex = data.loc[union,:]
del data_redex["sample"]
disease_match = {'kidney chromophobe':'KICH', 'colon adenocarcinoma':'COAD', 
                 'lung squamous cell carcinoma':'LUSC', 'bladder urothelial carcinoma':'BLCA', 
                 'diffuse large B-cell lymphoma':'DLBC', 'adrenocortical cancer':'ACC',
                 'cholangiocarcinoma':'CHOL', 'kidney papillary cell carcinoma':'KIRP', 
                 'esophageal carcinoma':'ESCA', 'liver hepatocellular carcinoma':'LIHC', 
                 'uterine carcinosarcoma':'UCS', 'glioblastoma multiforme':'GBM', 
                 'stomach adenocarcinoma':'STAD', 'uveal melanoma':'UVM', 
                 'pancreatic adenocarcinoma':'PAAD', 'thyroid carcinoma':'THCA', 
                 'sarcoma':'SARC', 'uterine corpus endometrioid carcinoma':'UCEC', 
                 'thymoma':'THYM', 'kidney clear cell carcinoma':'KIRC', 
                 'acute myeloid leukemia':'LAML', 'skin cutaneous melanoma':'SKCM', 
                 'lung adenocarcinoma':'LUAD', 'mesothelioma':'MESO', 
                 'brain lower grade glioma':'LGG', 'breast invasive carcinoma':'BRCA', 
                 'testicular germ cell tumor':'TGCT', 'rectum adenocarcinoma':'READ', 
                 'head & neck squamous cell carcinoma':'HNSC', 
                 'cervical & endocervical cancer':'CESC', 
                 'pheochromocytoma & paraganglioma':'PCPG', 
                 'prostate adenocarcinoma':'PRAD', 
                 'ovarian serous cystadenocarcinoma':'OV'}
primary_disease = []
for each in data_redex["_primary_disease"].tolist():
    primary_disease.append(disease_match[each])
data_redex["primary_disease"] = primary_disease
matrix_sample = tcga_matrix_t.index.tolist()
matrix_sample_correct = []
for each in matrix_sample:
    matrix_sample_correct.append(each.replace("-","."))
tcga_matrix_t.index = matrix_sample_correct
tcga_matrix_t = tcga_matrix_t.loc[union,:]
quan_propotion = quan_propotion.loc[union,:]
new_col = []
for each in quan_propotion.columns.tolist():
    new_col.append("QS_"+each)
quan_propotion.columns = new_col
LM_propotion = LM_propotion.loc[union,:]
LM_propotion = LM_propotion.dropna(axis=0,how='any')
new_col = []
for each in LM_propotion.columns.tolist():
    new_col.append("LM_"+each)
LM_propotion.columns = new_col
propotion = propotion.rename(columns = {"Unnamed: 0":"cellID"})
propotion.index = propotion["cellID"].tolist()
propotion = propotion.loc[union,:]
new_col = []
for each in propotion.columns.tolist():
    if each == "cellID":
        new_col.append("cellID")
    else:
        new_col.append("EPIC_"+each)
propotion.columns = new_col
columns_index_EPIC = propotion.columns.tolist()
columns_index_QS = quan_propotion.columns.tolist()
columns_index_LM = LM_propotion.columns.tolist()
tcga_obs = dict()
for each in columns_index_EPIC:
    tcga_obs[each] = propotion[each].tolist()
for each in columns_index_QS:
    tcga_obs[each] = quan_propotion[each].tolist()
for each in columns_index_LM:
    tcga_obs[each] = LM_propotion[each].tolist()
columns_index = data_redex.columns.tolist()
for each in columns_index:
    tcga_obs[each] = data_redex[each].tolist()
tcga_arr = np.array(tcga_matrix_t)
tcga_arr_t = np.array(tcga_matrix_t.T)
a = API.DatabaseAPI("tcga")
a.write_collection_matrix_obs_by_var(tcga_arr,overwrite=True)
a.write_collection_gene_matrix_var_by_obs(tcga_arr_t,genes,overwrite = True)
a.write_collection_obs(tcga_obs,overwrite=True)
tcga_var = dict()
tcga_var["geneSymbol"] = genes
a.write_collection_var(tcga_var,overwrite = True)

# gtex
gtex_matrix = pd.read_csv("./gtexTpm.csv")
genes = gtex_matrix["gene"].tolist()
gtex_matrix_t = gtex_matrix.iloc[:,1:].T
gtex_matrix_t.columns = genes
data = pd.read_csv("./GTEX_phenotype",sep="\t")
propotion = pd.read_csv("./gtex_EPIC_propotion.csv")
cellID = propotion["Unnamed: 0"].tolist()
quan_propotion = pd.read_csv("./gtex_QS_propotion.csv")
cell_type = quan_propotion["cell_type"].tolist()
quan_propotion = quan_propotion.iloc[:,2:].T
quan_propotion.columns = cell_type
LM_propotion = pd.read_csv("./gtex_LM_propotion.csv")
cols = LM_propotion["Unnamed: 0"].tolist()
cell_type = LM_propotion.columns.tolist()[1:]
LM_propotion = LM_propotion.iloc[:,1:]
LM_propotion.columns = cell_type
LM_propotion.index = cols
sample = data["Sample"].tolist()
sample_correct = []
for each in sample:
    sample_correct.append(each.replace("-","."))
data.index = sample_correct
union = [i for i in sample_correct if i in cellID]
data_redex = data.loc[union,:]
del data_redex["Sample"]
matrix_sample = gtex_matrix_t.index.tolist()
matrix_sample_correct = []
for each in matrix_sample:
    matrix_sample_correct.append(each.replace("-","."))
gtex_matrix_t.index = matrix_sample_correct
gtex_matrix_t = gtex_matrix_t.loc[union,:]
quan_propotion = quan_propotion.loc[union,:]
new_col = []
for each in quan_propotion.columns.tolist():
    new_col.append("QS_"+each)
quan_propotion.columns = new_col
LM_propotion = LM_propotion.loc[union,:]
new_col = []
for each in LM_propotion.columns.tolist():
    new_col.append("LM_"+each)
LM_propotion.columns = new_col
propotion = propotion.rename(columns = {"Unnamed: 0":"cellID"})
propotion.index = propotion["cellID"].tolist()
propotion = propotion.loc[union,:]
new_col = []
for each in propotion.columns.tolist():
    if each == "cellID":
        new_col.append("cellID")
    else:
        new_col.append("EPIC_"+each)
propotion.columns = new_col
columns_index_EPIC = propotion.columns.tolist()
columns_index_QS = quan_propotion.columns.tolist()
columns_index_LM = LM_propotion.columns.tolist()
gtex_obs = dict()
for each in columns_index_EPIC:
    gtex_obs[each] = propotion[each].tolist()
for each in columns_index_QS:
    gtex_obs[each] = quan_propotion[each].tolist()
for each in columns_index_LM:
    gtex_obs[each] = LM_propotion[each].tolist()
columns_index = data_redex.columns.tolist()
for each in columns_index:
    gtex_obs[each] = data_redex[each].tolist()
gtex_arr = np.array(gtex_matrix_t)
gtex_arr_t = np.array(gtex_matrix_t.T)
b = API.DatabaseAPI("gtex")
b.write_collection_matrix_obs_by_var(gtex_arr,overwrite=True)
b.write_collection_gene_matrix_var_by_obs(gtex_arr_t,genes,overwrite = True)
b.write_collection_obs(gtex_obs,overwrite=True)
gtex_var = dict()
gtex_var["geneSymbol"] = genes
b.write_collection_var(gtex_var,overwrite = True)

#EPIC
EPIC = pd.read_csv("./refTpm21765Genes.csv")
EPIC = EPIC.loc[EPIC["Unnamed: 0"].isin(genes),:]
genes_EPIC = EPIC["Unnamed: 0"].tolist()
EPIC_t = EPIC.iloc[:,1:].T
EPIC_t.columns = genes_EPIC
EPIC = EPIC_t.T
col_sum = EPIC.sum(axis = 0)
ax = 1e6 / col_sum
EPIC_TPM = EPIC * ax
EPIC_t = EPIC_TPM.T
EPIC_arr_t = np.array(EPIC_t.T)
genes = EPIC_t.columns.tolist()
EPIC_arr = np.array(EPIC_t)
b = API.DatabaseAPI("ref")
b.write_collection_matrix_obs_by_var(EPIC_arr,overwrite=True)
b.write_collection_gene_matrix_var_by_obs(EPIC_arr_t,genes,overwrite = True)

EPIC_var = dict()
EPIC_var["geneSymbol"] = genes
b.write_collection_var(EPIC_var,overwrite = True)
EPIC_obs = dict()
EPIC_obs["celltype"] = EPIC_t.index.tolist()
b.write_collection_obs(EPIC_obs,overwrite = True)

#LM22
LM22 = pd.read_csv("./LM22.txt", sep = "\t")
LM22 = LM22.loc[LM22["genesinput"].isin(genes),:]
genes_LM = LM22["genesinput"].tolist()
LM22_t = LM22.iloc[:,1:].T
LM22_t.columns = genes_LM
LM22 = LM22_t.T
col_sum = LM22.sum(axis = 0)
ax = 1e6 / col_sum
LM22_TPM = LM22 * ax
LM22_t = LM22_TPM.T
LM22_t.index = ['B.cells.naive',
 'B.cells.memory',
 'Plasma.cells',
 'T.cells.CD8',
 'T.cells.CD4.naive',
 'T.cells.CD4.memory.resting',
 'T.cells.CD4.memory.activated',
 'T.cells.follicular.helper',
 'T.cells.regulatory..Tregs.',
 'T.cells.gamma.delta',
 'NK.cells.resting',
 'NK.cells.activated',
 'Monocytes',
 'Macrophages.M0',
 'Macrophages.M1',
 'Macrophages.M2',
 'Dendritic.cells.resting',
 'Dendritic.cells.activated',
 'Mast.cells.resting',
 'Mast.cells.activated',
 'Eosinophils',
 'Neutrophils']
LM22_arr_t = np.array(LM22_t.T)
genes = LM22_t.columns.tolist()
LM22_arr = np.array(LM22_t)
b = API.DatabaseAPI("LM_ref")
b.write_collection_matrix_obs_by_var(LM22_arr,overwrite=True)
b.write_collection_gene_matrix_var_by_obs(LM22_arr_t,genes,overwrite = True)

LM22_var = dict()
LM22_var["geneSymbol"] = genes
b.write_collection_var(LM22_var,overwrite = True)
LM22_obs = dict()
LM22_obs["celltype"] = LM22_t.index.tolist()
b.write_collection_obs(LM22_obs,overwrite = True)

#QS
QS = pd.read_csv("./QS_TPM.csv")
genes_QS = QS["Unnamed: 0"].tolist()
QS_t = QS.iloc[:,1:].T
QS_t.columns = genes_QS
QS_t["label"] = QS.columns.tolist()[1:]
QS_t["label"] = ['B cells', 'B cells', 'B cells', 'B cells', 'B cells', 'B cells', 
                 'Macrophages M1', 'Macrophages M1', 'Macrophages M1', 'Macrophages M2', 
                 'Macrophages M2', 'Macrophages M2', 'Monocytes', 'Monocytes', 'Monocytes', 
                 'Monocytes', 'Monocytes', 'Monocytes', 'Neutrophils', 'Neutrophils', 'Neutrophils', 
                 'Neutrophils', 'Neutrophils', 'Neutrophils', 'NK cells', 'NK cells', 'NK cells', 
                 'NK cells', 'NK cells', 'NK cells', 'T cells CD4', 'T cells CD4', 'T cells CD4', 
                 'T cells CD4', 'T cells CD8', 'T cells CD8', 'T cells CD8', 'T cells CD8', 
                 'Tregs', 'Tregs', 'Tregs', 'Tregs', 'Tregs', 'Tumor', 'Tumor', 'Dendritic cells', 
                 'Dendritic cells', 'Dendritic cells', 'Dendritic cells', 'Dendritic cells', 
                 'Dendritic cells', 'Dendritic cells', 'Dendritic cells']
QS_t = QS_t.groupby(by = "label").mean()
QS_TPM = QS_t.T
QS_TPM = QS_TPM.loc[QS_TPM.index.isin(genes),:]
col_sum = QS_TPM.sum(axis = 0)
ax = 1e6 / col_sum
QS_TPM = QS_TPM * ax
QS_t = QS_TPM.T
QS_t.index = ['B cell', 'Myeloid dendritic cell', 'Macrophage M1', 'Macrophage M2', 'Monocyte',
       'NK cell', 'Neutrophil', 'T cell CD4+(non-regulatory)', 'T cell CD8+', 
    'T cell regulatory(Tregs)', 'Tumor']
QS_arr_t = np.array(QS_t.T)
genes = QS_t.columns.tolist()
QS_arr = np.array(QS_t)
b = API.DatabaseAPI("QS_ref")
b.write_collection_matrix_obs_by_var(QS_arr,overwrite=True)
b.write_collection_gene_matrix_var_by_obs(QS_arr_t,genes,overwrite = True)
QS_var = dict()
QS_var["geneSymbol"] = genes
b.write_collection_var(QS_var,overwrite = True)
QS_obs = dict()
QS_obs["celltype"] = QS_t.index.tolist()
b.write_collection_obs(QS_obs,overwrite = True)
