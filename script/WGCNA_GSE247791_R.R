
# RNA-seq WGCNA Analysis (GSE247791)
# Najneen Rejwana


setwd("F:/github/WGCNA/WCGNA_GSE247791")


## LOAD LIBRARIES---------------------------------------------------------------
library(DESeq2)
library(WGCNA)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()


## 1. LOAD DATA-----------------------------------------------------------------
count_data <- read.csv("count_data_GSE247791_EDIT.csv", row.names = 1)
meta_data  <- read.csv("metadata_EDIT_GSE247791_csv.csv", row.names = 1)

# Align samples
count_data <- count_data[, rownames(meta_data)]
stopifnot(all(colnames(count_data) == rownames(meta_data)))


## 2. DESEQ2 NORMALIZATION (VST)------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = meta_data,
                              design = ~ treatment)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# Samples as rows for WGCNA
datExpr0 <- t(assay(vsd))


## 3. DATA CLEANING (GOOD GENES/SAMPLES)----------------------------------------
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr0[gsg$goodSamples, gsg$goodGenes]
} else {
  datExpr <- datExpr0
}


## 4. SAMPLE CLUSTERING (OUTLIERS)----------------------------------------------
jpeg("Sample_clustering.jpg", width=800, height=600)
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering", xlab = "", sub = "")
dev.off()


## 5. SOFT-THRESHOLD PLOT (MANUAL SELECTION)------------------------------------
powers <- c(1:10, seq(12,20,2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

jpeg("Soft_threshold.jpg", width=1000, height=500)
par(mfrow=c(1,2))
# Scale-free topology
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit (R²)",
     type="n")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.8, col="blue")
# Mean connectivity
plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n")
text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels=powers, col="red")
dev.off()

#### MANUALLY change this based on plot:
softPower <- 14   # <-- set according to R² ≥ 0.8
##Here, check the figure named Soft_threshold.jpg. softPower <- 14 can be changed based on figure generated. Here, Scale-Free Fit: At power 12, your $R^2$ is just below the 0.8 line (looks like ~0.76). At power 14, the value officially crosses or sits right on that 0.8 threshold.


## 6. NETWORK CONSTRUCTION------------------------------------------------------
cor <- WGCNA::cor  # fix any conflict
net <- blockwiseModules(datExpr,
                        power = softPower,
                        TOMType = "unsigned",
                        minModuleSize = 30,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        verbose = 3)

moduleColors <- labels2colors(net$colors)

## 7. MODULE DENDROGRAM---------------------------------------------------------
jpeg("Gene_dendrogram_modules.jpg", width=900, height=700)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03)
dev.off()


## 8. MODULE–TRAIT RELATIONSHIP--------------------------------------------------
traitData <- meta_data
traitData$treatment_num <- ifelse(traitData$treatment=="Abemaciclib",1,0)
traitData <- traitData[rownames(datExpr), , drop=FALSE]

MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)

moduleTraitCor <- cor(MEs, traitData$treatment_num, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Increase height slightly to accommodate the bottom labels
jpeg("Module_trait_heatmap.jpg", width=900, height=1600, res=120)

# Adjust margins: c(bottom, left, top, right)
par(mar = c(10, 15, 3, 3)) 

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "Abemaciclib Treatment",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = signif(moduleTraitCor, 2),
               cex.text = 1.0,
               zlim = c(-1,1),
               main = "Module–Trait Relationship (GSE247791)",
               
               # --- FIXING X-AXIS ALIGNMENT ---
               xLabelsAngle = 0,      # Keeps the text horizontal (0 degrees)
               xLabelsAdj = 0.5,      # Centers the text under the column
               # -------------------------------
               
               setStdMargins = FALSE)

dev.off()

## 9. SELECT BEST MODULE--------------------------------------------------------

moduleTrait <- data.frame(
  module = names(MEs),
  cor = as.numeric(moduleTraitCor),
  p = as.numeric(moduleTraitPvalue)
)
moduleTrait <- moduleTrait[order(moduleTrait$p), ]
print(moduleTrait)

bestModule <- gsub("ME","",moduleTrait$module[1])
print(paste("Best module:", bestModule))


write.csv(moduleTrait, 
          file = "GSE247791_Module_Trait_Summary.csv", 
          row.names = FALSE)


## 10. EXTRACT GENES & SAVE-----------------------------------------------------------------
moduleGenes <- colnames(datExpr)[moduleColors==bestModule]
write.csv(moduleGenes, paste0("genes_", bestModule, "_module.csv"), row.names=FALSE)

geneModule <- data.frame(Gene=colnames(datExpr), Module=moduleColors)
write.csv(geneModule, "GSE247791_WGCNA_gene_modules.csv", row.names=FALSE)

