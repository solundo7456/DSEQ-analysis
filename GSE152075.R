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
meta <- pData(gse)
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
unique(meta2$age)
meta2$age[meta2$age %in% c("NA")] <- NA 
meta2$age <- as.numeric(meta2$age)

unique(meta2$viral_load)
meta2$viral_load[meta2$viral_load %in% c("NA")] <- NA 
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

#1. functionality
# Install necessary packages if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)

# 2. Filter significant genes
# -----------------------------
sig_genes <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)]
length(sig_genes)
head(sig_genes)
sig_genes
# 3. Map gene symbols to Entrez IDs
# -----------------------------
entrez_ids <- bitr(sig_genes,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
head(entrez_ids)

# 4. GO enrichment (Biological Process)
# -----------------------------
ego <- enrichGO(gene         = entrez_ids$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod= "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)

ego_1 <- enrichGO(gene         = entrez_ids$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "MF",
                pAdjustMethod= "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)
# Top GO terms
head(ego)
head(ego_1)
# Dotplot
dotplot(ego, showCategory=20) + ggtitle("GO Biological Process Enrichment")

# 5. KEGG pathway enrichment
# -----------------------------
options(timeout = 600) 
ekegg <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                    organism     = "hsa",
                    pvalueCutoff = 0.05)

# Top KEGG pathways
head(ekegg)
dotplot(ekegg, showCategory=20) + ggtitle("KEGG Pathway Enrichment")

# -----------------------------
# 6. Reactome pathway enrichment
# -----------------------------
ereact <- enrichPathway(gene=entrez_ids$ENTREZID,
                        organism="human",
                        pvalueCutoff=0.05,
                        readable=TRUE)

# Top Reactome pathways
head(ereact)
dotplot(ereact, showCategory=20) + ggtitle("Reactome Pathway Enrichment")
