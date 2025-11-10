#creating our table

counts <- read.table("GSE152075_raw_counts_GEO.txt",
		     header = TRUE,
		    row.names = 1,
		    sep = ' ',
		   check.names = FALSE)
head(counts)

library(GEOquery)
library(limma)
library(statmod)

#getting our metadata

file_path <- "GSE152075_series_matrix.txt.gz"
# Load the GEO series matrix
gse <- getGEO(filename = file_path)
#understanding the metadata
class(gse)
head(gse)
str(gse)
slotNames(gse)
#inspecting Expression values 
exprs(gse)
#sample metadata
head(pData(gse))
#feature metadata
fData(gse)
#experimental metadata
experimentData(gse)

head(counts)
dim(counts)
#create sample list
meta <- pData(gse[[1]])
head(meta)
#ensure columns of counts match rowa of metadata
all(colnames(counts) %in% rownames(meta$title))
#confirm which is not present
setdiff(colnames(counts),rownames(meta$title))


#library sizes
library_sizes <- colSums(counts)
summary(library_sizes)
hist(log10(library_sizes), main='library sizes (log10)', xlab='log10(Total counts)')

#filter remove genes with low expression
keep <- rowSums(counts >10) >=5
counts_filtered <- counts[keep,]
dim(counts_filtered)

#prepare metadata 
#reorder
meta2 <- meta[match(colnames(counts_filtered), meta$title), ]
all(meta2$title == colnames(counts_filtered))
meta2
str(meta2)
#renaming 
colnames(meta2) <- gsub(":ch1", "", colnames(meta2))
meta2$infection_status <- meta2$'sars-cov-2 positivity'
meta2$viral_load <- meta2$n1_ct
str(meta2)

#data type
meta2$infection_status <- factor(meta2$infection_status)
meta2$`sequencing_batch` <- factor(meta2$sequencing_batch)
meta2$gender <- factor(meta2$gender)
#handling numeric
meta2$age <- as.numeric(meta2$age)
meta2[is.na(meta2$age), "age"]
meta2$age[meta2$age %in% c("NA")] <- NA 
meta2$age <- as.numeric(meta2$age)
meta2$viral_load <- as.numeric(meta2$viral_load)

#reorder
head(colnames(counts_filtered))
head(rownames(meta2))
rownames(meta2) <- meta2$title 
meta2 <- meta2[colnames(counts_filtered), ]
#creating DEseq2

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = meta2, design = ~infection_status)
class(dds)
#run analysis
dds <- DESeq(dds)
plotDispEsts(dds)
#results
res <- results(dds, alpha=0.05)
summary(res)

#number of significate genes
sum(res$padj < 0.05, na.rm = TRUE)  # number of DE genes at FDR 5%

#MA-plot
plotMA(res, main="DESeq2 MA-plot")

#volcano plot
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Positive vs Negative")
#Heatmap 
library(pheatmap)
top_genes <- head(order(res$padj), 50)
# top 50 DE genes
mat <- assay(vst(dds))[top_genes, ]
colnames(meta2)

pheatmap(mat, scale="row", annotation_col = meta2[, c("infection_status","gender")])
write.csv(as.data.frame(res), "DE_results_Positive_vs_Negative.csv")

