# CIBERSORT R script v1.04 (last updated 10-24-2016)
# Note: Signature matrix construction is not currently available; use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('e1071')
#       install.pacakges('parallel')
#       install.packages('preprocessCore')
#       if preprocessCore is not available in the repositories you have selected, run the following:
#           source("http://bioconductor.org/biocLite.R")
#           biocLite("preprocessCore")
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
# This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('CIBERSORT.R')
#       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN, absolute, abs_method)
#
#       Options:
#       i)   perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#       ii)  QN = Quantile normalization of input mixture (default = TRUE)
#       iii) absolute = Run CIBERSORT in absolute mode (default = FALSE)
#               - note that cell subsets will be scaled by their absolute levels and will not be
#                 represented as fractions (to derive the default output, normalize absolute
#                 levels such that they sum to 1 for each mixture sample)
#               - the sum of all cell subsets in each mixture sample will be added to the ouput
#                 ('Absolute score'). If LM22 is used, this score will capture total immune content.
#       iv)  abs_method = if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
#               - sig.score = for each mixture sample, define S as the median expression
#                 level of all genes in the signature matrix divided by the median expression
#                 level of all genes in the mixture. Multiple cell subset fractions by S.
#               - no.sumto1 = remove sum to 1 constraint
#
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt


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
