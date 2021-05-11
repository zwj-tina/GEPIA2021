library(EPIC)
#gtex
bulkSamplesMatrix <- read.csv("./gtexTpm.csv")
rownames(bulkSamplesMatrix) <- bulkSamplesMatrix$gene
bulkSamplesMatrix <- bulkSamplesMatrix[,-1]
gtex <- EPIC(bulk = bulkSamplesMatrix)
#tcga
bulkSamplesMatrix <- read.csv("./tcgaTpm.csv")
rownames(bulkSamplesMatrix) <- bulkSamplesMatrix$gene
bulkSamplesMatrix <- bulkSamplesMatrix[,-1]
gtex <- EPIC(bulk = bulkSamplesMatrix)
