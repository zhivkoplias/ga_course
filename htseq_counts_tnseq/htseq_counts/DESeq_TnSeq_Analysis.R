source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)

directory <- "/home/erik/sweden/courses/2nd_semester/Genome_Analysis/ga_course/htseq_counts_tnseq/htseq_counts/"

sampleFiles_Serum <- grep("Serum",list.files(directory),value=TRUE)
sampleCondition_Serum <- sub("(.*Serum).*","\\1",sampleFiles_Serum)
sampleTable_Serum <- data.frame(sampleName = sampleFiles_Serum,
                          fileName = sampleFiles_Serum,
                          condition = sampleCondition_Serum)

sampleFiles_Human <- grep("Human",list.files(directory),value=TRUE)
sampleCondition_Human <- sub("(.*Human).*","\\1",sampleFiles_Human)
sampleTable_Human <- data.frame(sampleName = sampleFiles_Human,
                                fileName = sampleFiles_Human,
                                condition = sampleCondition_Human)

sampleFiles_BH <- grep("BH",list.files(directory),value=TRUE)
sampleCondition_BH <- sub("(.*BH).*","\\1",sampleFiles_BH)
sampleTable_BH <- data.frame(sampleName = sampleFiles_BH,
                                fileName = sampleFiles_BH,
                                condition = sampleCondition_BH)

sampleCondition <- c(sampleCondition_Serum, sampleCondition_BH, sampleCondition_Human)

sampleTable_SH <- merge(sampleTable_BH, sampleTable_Serum, all=TRUE)
sampleTable <- merge(sampleTable_SH, sampleTable_Human, all=TRUE)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

#prefiltering

keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]
dds <- ddsHTSeq 

#Note on factor levels

#dds$condition <- factor(dds$condition, levels = c("Serum","Human"))
#dds$condition <- relevel(dds$condition, ref = "Serum")
#dds$condition <- droplevels(dds$condition)


# DE analysis
dds <- dds[, dds$condition %in% c("Serum","BH","Human") ]
dds <- DESeq(dds)
res <- results(dds)
res

res <- results(dds, name="condition_Serum_vs_BH")


res <- results(dds, contrast=c("condition","Serum","BH"))

# Log fold change shrinkage for visualization and ranking
#biocLite("apeglm")
library(apeglm)

resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_Serum_vs_BH", type="apeglm")
resLFC

# p-values ordering

resOrdered <- res[order(res$log2FoldChange),]
summary(res)

sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)

#plots DE

plotMA(res, ylim=c(-9,9))

plotMA(resLFC, ylim=c(-9,9))

idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

# plots counts

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.ontop = TRUE
  )

# pca

# install.packages("factoextra")
library(factoextra)

setwd("/home/erik/sweden/courses/2nd_semester/Genome_Analysis/ga_course/htseq_counts_tnseq/htseq_counts/for_pca")
temp = list.files(pattern="*.txt")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i], header=TRUE, sep="\t"))

BH1 <- merge(BHI_ERR1801012_transcripts.txt, BHI_ERR1801013_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
BH2 <- merge(BH1, BHI_ERR1801014_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
H1 <- merge(BH2, Human_ERR1801009_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
H2 <- merge(H1, Human_ERR1801010_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
H3 <- merge(H2, Human_ERR1801011_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
S1 <- merge(H3, Serum_ERR1801006_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
S2 <- merge(S1, Serum_ERR1801007_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
S3 <- merge(S2, Serum_ERR1801008_transcripts.txt, by.x="aaaGene", by.y="aaaGene")

dataset_for_pca <- S3
dataset_for_pca[, 2:10] <- sapply(dataset_for_pca[, 2:10], as.numeric)
dataset_for_pca <- dataset_for_pca[, 2:10]
str(dataset_for_pca)
#install.packages("ggfortify")

res.pca = prcomp(t(as.matrix(dataset_for_pca)))
project.pca <- res.pca
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
# attemot 1
project.pca$x
autoplot(project.pca)

require(scatterplot3d)
par(mar=c(4,4,4,4), cex=1.0, cex.main=0.8, cex.axis=0.8)

scatterplot3d(project.pca$x[,1:3], angle=-40, main="", color="black", pch=17, xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"), zlab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"), grid=FALSE, box=FALSE)
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(project.pca$x[,1:3], grid = c("xy", "xz", "yz"))

source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(project.pca$x[,1:3], grid = c("xy", "xz", "yz"))

# PCA with colour

library(FactoMineR)
iris.pca <- PCA(iris, quali.sup=5)
plot(project.pca, habillage = 5, 
     col.hab = c("green", "blue", "red"), 
     title = "Dataset projected onto PC1-2 Subspace")

#Another attempt

category = c( "BHI","BHI","BHI","Human","Human","Human","Serum","Serum","Serum")
qplot(project.pca$x[,1],project.pca$x[,2], colour=category, size=9, xlab="PCA 1", ylab="PCA 2") 

#Third
df_out <- as.data.frame(project.pca$x)
df_out$group <- sapply( strsplit(as.character(colnames(dataset_for_pca)), "_"), "[[", 1 )
head(df_out)

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point()
p

# Claudios - on DE data

rld <- rlogTransformation(dds, blind=TRUE)
plotPCA(rld)

# more results

mcols(res)$description


#save data

save_it <- as.data.frame(resOrdered)

write.csv(as.data.frame(save_it[1]), 
          file="exclude_two_sample_3_cond_down.csv")

new <- read.csv("exclude_two_sample_3_cond_down.csv", sep = ',')
names <- new$X
write.csv(names, 
          file="exclude_two_sample_3_cond_gene_names_down.csv", row.names = FALSE)
