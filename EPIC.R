library(EPIC)
bulkSamplesMatrix <- read.csv("../gtexTpmGenes.csv")
rownames(bulkSamplesMatrix) <- bulkSamplesMatrix$gene
bulkSamplesMatrix <- bulkSamplesMatrix[,-1]
gtex <- EPIC(bulk = bulkSamplesMatrix)

bulkSamplesMatrix <- read.csv("../tcgaTpmGenes.csv")
rownames(bulkSamplesMatrix) <- bulkSamplesMatrix$gene
bulkSamplesMatrix <- bulkSamplesMatrix[,-1]
gtex <- EPIC(bulk = bulkSamplesMatrix)
