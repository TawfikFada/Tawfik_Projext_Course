# Name: Tawfik Fada 
# LiU ID: Tawfa921
# Email: Tawfa921@student.liu.se
# Last edited: 17/12/2023
# Script description: Visualization and differential expression analysis of Mock and Treated samples. 
#-------------------------------------------------------------------------------
# Setting working directory
setwd("~/Prjoect course R Data")
#-------------------------------------------------------------------------------
# Loading packages 
library(dplyr)
library(ggplot2)
library(tidyverse)
#-------------------------------------------------------------------------------
file.path("Galaxy166-[edgeR_normcounts].tabular")
file.path("Galaxy168-[Galaxy168-[annotateMYIDs_on_data_165__Annotated_IDs].tabular")
#-------------------------------------------------------------------------------
# Reading the annotated data derived from galaxy
annotation <- read.table("annotation.tabular", header = TRUE)

# Reading normcounts data
gene_expression <- read.table("edgeR_normcounts.tabular", header = TRUE)
#-------------------------------------------------------------------------------
# Creating table to show all of our expressed genes
gene_expression <- merge(gene_expression, annotation, by = 1, all.x = TRUE)
#-------------------------------------------------------------------------------
# Looking at the data to learn what I am working with a bit better
dim(gene_expression)
head(gene_expression)
tail(gene_expression)
summary(gene_expression)
#-------------------------------------------------------------------------------
# Removing NA from gene_expression
table(is.na(gene_expression))
gene_expression <- na.omit(gene_expression)
#-------------------------------------------------------------------------------
# Changing row names to symbol 
rownames(gene_expression) <- gene_expression$SYMBOL

gene_expression$GeneID <- NULL
gene_expression$SYMBOL <- NULL
#-------------------------------------------------------------------------------
# Sorting the data in alphabetical order
gene_expression <- gene_expression[order(rownames(gene_expression)),]
#-------------------------------------------------------------------------------
# Creating the meta data files
meta_data <- data_frame(colnames(gene_expression))
colnames(meta_data) <- "SampleID"
meta_data$group <- rep(c("treated","mock"), times=c(3,3))
meta_data$sampletype <- "A595"
#-------------------------------------------------------------------------------
library(tidyverse)
# calculating principle components
count_pca <- prcomp(t(gene_expression))
#-------------------------------------------------------------------------------
library(ggfortify)# needed to run autoplot
# Visualizing results
autoplot(count_pca)

autoplot(count_pca, data=meta_data, colour="group", size=4)

autoplot(count_pca, data=meta_data, colour="group", size=4, frame=TRUE, frame.type = 'norm')+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#-------------------------------------------------------------------------------
# Performing a T-test
pvalue <- data.frame(p_value=rep(0,nrow(gene_expression)))
# Renaming the rows of pvalue to match that of gene_expression
rownames(pvalue) <- rownames(gene_expression)
# Using a loop function to insert the p-value for each row into the table of p-values.
for(i in 1:nrow(gene_expression)){
  pvalue$p_value[i]<- t.test(gene_expression [i,1:3], gene_expression[i,4:6],
                             alternative = "two.sided")$p.value 
}
# adjusting the P-values
pvalue$padjust <- p.adjust(pvalue$p_value, "BH")
# logFC (Fold Change) is employed for additional gene filtering, given the initial count of 4238 genes is too extensive.
pvalue$logFC <- (rowMeans(gene_expression[, 1:3]))-rowMeans((gene_expression[, 4:6]))
# Filtering values below 0.05
significant_pvalue <- pvalue[pvalue$padjust<0.05,,drop=FALSE]

significant_logFC <- significant_pvalue[significant_pvalue$logFC > 2 | significant_pvalue$logFC < -2,,drop=FALSE]
#-------------------------------------------------------------------------------
# Plot #1: Volcano Plot

# Basic volcano plot
p <- pvalue %>%
  ggplot(mapping = aes(x=logFC, y=-log10(padjust)))+
  geom_point()+
  theme_minimal()

p

# Incorporate a line to emphasize the padjusted value of 0.05 and establish a threshold for the logFC.
p2 <- p +
  geom_hline(yintercept = -log10(0.05), col="red") +
  geom_vline(xintercept= c(-2, 2), col="blue")

p2

#Assign colors to the dots according to their differential expression.
pvalue$diffexpressed <- "NO"
pvalue$diffexpressed[pvalue$logFC > 2 & pvalue$padjust < 0.05] <- "UP"
pvalue$diffexpressed[pvalue$logFC < -2 & pvalue$padjust < 0.05] <- "DOWN"

p <- pvalue %>%
  ggplot(mapping=aes(x=logFC, y=-log10(padjust), col=diffexpressed))+
  geom_point()+
  theme_minimal()

p2 <- p +
  geom_hline(yintercept = -log10(0.05), col="red")+
  geom_vline(xintercept= c(-2,2), col="blue")

mycolors <- c("orange", "purple", "darkgrey")

names(mycolors) <- c("DOWN", "UP", "NO") 

p3 <- p2+
  scale_color_manual(values=mycolors)

p3
#-------------------------------------------------------------------------------
# Plot #2: Scatter plot
library(ggrepel)
library(pheatmap)

# Generate a data frame named FDR_pvalue that includes adjusted p-values and logFC for all 13,000 genes, 
# with gene ids as row names.
FDR_pvalue <- data.frame(pvalue$padjust)
rownames(FDR_pvalue) <- rownames(gene_expression)

FDR_pvalue

# Creating significant_pvalue with keep FDR p value<0.05 
significant_pvalue2 <- FDR_pvalue[FDR_pvalue$padjust<0.05,,drop=FALSE]

# adding the logFC values into the FDR_pvalue data frame
logFC <- (rowMeans(gene_expression[i, 1:3])-rowMeans(gene_expression[i, 4:6]))

FDR_pvalue$logFC <- 0

FDR_pvalue$logFC <- logFC 

# Generate a data frame named top_genes, containing the top 7 up-regulated and top 7 down-regulated genes.

top_genes <- significant_pvalue[c(1:7, 4231:4238),]

top_genes <- rownames(top_genes) 

top_genes_expression <- gene_expression[(rownames(gene_expression) %in% top_genes),]


top_genes_expression$MockMean<-rowMeans(top_genes_expression[,1:3])  
top_genes_expression$TreatedMean<-rowMeans(top_genes_expression[,4:6]) 

head(top_genes_expression)

top_genes_expression$genename<-rownames(top_genes_expression)

# Scatter plot with the top 7 Up and Down regulated genes 
ggplot(top_genes_expression, aes(x=MockMean, y=TreatedMean, label=genename))+
  geom_point()+
  geom_abline(intercept = 0)+
  theme_bw()+
  ggtitle(label="Top 7 Up and Down Regulated Genes")+
  aes(color=TreatedMean)+
  geom_text_repel()+
  theme(plot.title = element_text(hjust = 0.5, size=20))
#-------------------------------------------------------------------------------
# Plot 3: Heatmap
pheatmap(top_genes_expression[1:8])
#-------------------------------------------------------------------------------
# Pathway analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# Generate a universe file that includes all the measured genes.
universe <- data.frame(gene=rownames(gene_expression)) 
universe <- bitr(universe$gene, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db" ) 

table(duplicated(universe$ENTREZID)) 
universe <- data.frame(gene=unique(universe$ENTREZID)) 

#Create a dataframe filtering for FDR values less than 0.05 and logFC values greater than 2 or less than -2.
differential_genes <- significant_logFC[significant_logFC == "NO"] <- NA 

differential_genes <- na.omit(significant_logFC)

# Define input parameters
significant_genes_pathway <- rownames(differential_genes)

genes_entrez <- bitr(significant_genes_pathway, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db) 

table(duplicated(genes_entrez$ENTREZID))

# KEGG pathway analysis
kegg1 <- enrichKEGG(gene= genes_entrez$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg",
                    universe = universe$ENTREZID, 
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH")


kegg1 <- setReadable(kegg1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(kegg1, showCategory=10, font.size=15) 

summary_enrich <- kegg1@result

# GO pathway analysis
go1 <- enrichGO(gene= genes_entrez$ENTREZID, 
                OrgDb = 'org.Hs.eg.db', 
                ont = "BP", 
                pvalueCutoff = 0.05, 
                pAdjustMethod = "BH", 
                universe = universe$ENTREZID, 
                qvalueCutoff = 0.2, 
                minGSSize = 10, 
                maxGSSize = 500, readable = TRUE)

dotplot(go1, showCategory=10, font.size=15) 

summary_go <- go1@result
