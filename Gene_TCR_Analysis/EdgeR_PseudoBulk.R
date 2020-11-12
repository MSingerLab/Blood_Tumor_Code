library(edgeR)
library(statmod)

# Load the matrix as a csv file
mtx <- 'Outputs/MousePseudoBulk_forEdgeR.csv'
mtx <- read.csv(mtx, sep = ",", row.names = 1, stringsAsFactors = FALSE)

# Convert the data.frame into a matrix
mtx <- as.matrix(t(mtx))

# Create a DGEList object
dgList <- DGEList(counts=mtx, genes = row.names(mtx))

# Filtering genes
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1

keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]
summary(cpm(dgList))

# Normalization
dgList <- calcNormFactors(dgList, method="TMM")

# Set up design matrix
Clone <- as.character(rep(1:110, 2))
Tissue <- rep("Blood", ncol(dgList))
Tissue[grep("Tumor", colnames(dgList))] <- "Tumor"

designMat <- model.matrix(~Clone+Tissue)
row.names(designMat) <- colnames(dgList)

# Estimating Dispersion
dgList <- estimateDisp(dgList, designMat, robust = TRUE)
dgList$common.dispersion

# DE genes
fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit)
topTags(lrt)

# Checks
o <- order(lrt$table$PValue)
cpm(dgList)[o[1:10],]

summary(decideTests(lrt))

plotMD(lrt)
abline(h=c(-1, 1), col="blue")

# Save results
edgeR_result <- topTags(lrt, n=length(row.names(lrt$table)))
DEgenes <- edgeR_result$table[edgeR_result$table$FDR < 0.05,]
DEgenes

# The output is 'unsorted' because in the Python script, we will
# sort the output before saving the final DE gene lists
write.csv(x=DEgenes, file='Outputs/MousePseudoBulk_DEgenes_unsorted.csv', row.names = F)

