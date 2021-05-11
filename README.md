# GEPIA2021
## process
    pre-process.py
        Preprocess data of Tcga and Gtex, transform to TPM format
    CIBERSORT_process.R
        CIBERSORT is an analytical tool from the Alizadeh Lab developed by Newman et al. to provide an estimation of the abundances of member cell types in a mixed cell        population, using gene expression data.
        This file is used to calculate the abundances of member cell types in a mixed cell population of Tcga and Gtex datasets by CIBESORT.
    EPIC.R
        Package implementing EPIC method to estimate the proportion of immune, stromal, endothelial and cancer or other cells from bulk gene expression data. It is based on reference gene expression profiles for the main non-malignant cell types and it predicts the proportion of these cells and of the remaining "other cells" (that are mostly cancer cells) for which no reference profile is given.
    This file is used to predicts the proportion of cells of Tcga and Gtex datasets.
    quantiseq.R
        quanTIseq is a computational pipeline for the quantification of the Tumor Immune contexture from human RNA-seq data.
        This file is used to calculate the quantification of the Tumor Immune contexture of Tcga and Gtex datasets.
## result
    tcga_QS_propotion.csv : calculated by quantiseq 
    tcga_LM_propotion.csv : calculated by CIBERSORT
    tcga_EPIC_propotion.csv : calculated by EPIC
    gtex_QS_propotion.csv : calculated by quantiseq
    gtex_LM_propotion.csv : calculated by CIBERSORT
    gtex_EPIC_propotion.csv : calculated by EPIC
## backend
    API.py
        Read data from MongoDB database or write data to MondoDB.
    views.py
        Returns the corresponding value for the front end, including propotion, expression and survival.
    database.py
        Write results data into database.
## reference
    HUGO.txt : protein-coding genes.
    LM22.txt : reference dataset of CIBERSORT.
