source('CIBERSORT.R')
sig_matrix <- read.csv("./LM22.txt", sep = "\t")
rownames(sig_matrix)  <- sig_matrix$Gene.symbol
sig_matrix <- sig_matrix[,-1]

mixture_file2 <- read.csv("./gtexTpmGenes.csv")
rownames(mixture_file2) <- mixture_file2$gene
mixture_file2 <- mixture_file2[,-1]
LM_GTEX_part <- CIBERSORT(sig_matrix = sig_matrix, mixture_file = mixture_file2, perm=0, QN=FALSE, absolute=TRUE, abs_method='sig.score')
write.csv(LM_GTEX, "./gtex_LM_propotion.csv")

mixture_file <- read.csv("./tcgaTpmGenes.csv")
rownames(mixture_file) <- mixture_file$gene
mixture_file <- mixture_file[,-1]
LM_TCGA <- CIBERSORT(sig_matrix = sig_matrix, mixture_file = mixture_file, perm=0, QN=FALSE, absolute=TRUE, abs_method='sig.score')
write.csv(LM_TCGA, "./tcga_LM_propotion.csv")
