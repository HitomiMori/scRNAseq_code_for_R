# scRNAseq_code_for_R
Single-cell analysis of ER+ breast cancer PDX models

#Data and package loading

setwd("~/Documents/COH_2018-2020/SC31/SC31_scRNAseq_2020")

setwd("~/Documents/COH_2018-2020/GS3_scRNA_2019_2020/GS3_intagration")

library(Seurat)

library(dplyr)

library(ggplot2)

#SC31

data <- readRDS("scRNA_human_phase0.3.rds")

#GS3 

data <- readRDS("GS3_integrated.rds")

#
PDX <- data
remove(data)

#without_integration

DefaultAssay(PDX) <- "RNA"

#integration

DefaultAssay(PDX) <- "integrated"

#reset

Idents(PDX) <- "seurat_clusters"


#Number of cells including the data set

stat <- as.data.frame(table(PDX@active.ident, PDX@meta.data$Phase, PDX@meta.data$orig.ident))
colnames(stat) <- c("Cluster", "CellCycle", "Treatment", "CellNumber")
write.csv(stat, "072219.csv")
write.csv(stat, "SC31_05032021.csv")


# Analysis_1
UMAPPlot(PDX, pt.size=0.1)+ggtitle("Clusters")

UMAPPlot(PDX, group.by="orig.ident",cols=c("#0066FF", "#FF99FF"), pt.size=0.1)+ggtitle("Treatment")

UMAPPlot(PDX, group.by="Phase", pt.size=0.1)+ggtitle("Cell cycle")

UMAPPlot(PDX, pt.size=0.1, split.by="orig.ident")+ggtitle("Clusters")

UMAPPlot(PDX, group.by="orig.ident",cols=c("#0066FF", "#FF99FF"), pt.size=0.1, split.by="orig.ident")+ggtitle("Treatment")

UMAPPlot(PDX, group.by="Phase", pt.size=0.1, split.by="orig.ident")+ggtitle("Cell cycle")

# Analysis_2
#Cell_Number

stat <- as.data.frame(table(PDX@active.ident, PDX@meta.data$Phase, PDX@meta.data$orig.ident))
colnames(stat) <- c("Cluster", "CellCycle", "Treatment", "CellNumber")
write.csv(stat, "SC31_06082021_8Cluster.csv")

#bar_Graph

stat <- as.data.frame(table(PDX@meta.data$Phase, PDX@meta.data$orig.ident))
colnames(stat) <- c("Phase", "Treatment", "Frequency")

g1 <- ggplot(data=stat, aes(x=Treatment, y=Frequency, fill=Phase)) +
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()
g1

stat <- as.data.frame(table(PDX@active.ident, PDX@meta.data$orig.ident))
colnames(stat) <- c("Cluster", "Treatment", "Frequency")

g2 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Treatment)) +
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()
g2

stat <- as.data.frame(table(PDX@active.ident, PDX@meta.data$Phase))
colnames(stat) <- c("Cluster", "Phase", "Frequency")

g3 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Phase)) +
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()
g3

stat <- as.data.frame(table(PDX@active.ident, PDX@meta.data$orig.ident, PDX@meta.data$Phase))
colnames(stat) <- c("Cluster", "Treatment", "Phase", "Frequency")

g4 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Phase)) +
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  facet_grid( ~ Treatment) 
g4

#
Gene <- "IL24"

sample.info <- data.frame(row.names = attr(PDX@active.ident, 'names'))
indexGene <- match(Gene, rownames(PDX@assays$RNA@counts))
sample.info$GENE = 'Neg'
sample.info$GENE[which(PDX@assays$RNA@counts[rownames(PDX@assays$RNA@counts)[indexGene], ] != 0)] <- 'Pos'
PDX <- AddMetaData(object = PDX, metadata = sample.info)
UMAPPlot(PDX, group.by="GENE", cols=c("lightgrey", "tomato"), pt.size = 0.1)+
  ggtitle(Gene) 
#
Idents(PDX) <- "GENE"
sub <- subset(PDX, idents = "Pos") 
UMAPPlot(sub,cols = "tomato")+ggtitle("IL24+")
Idents(sub) <- "seurat_clusters"
UMAPPlot(sub)+ggtitle("IL24+_Clusters")

stat <- as.data.frame(table(sub@active.ident, sub@meta.data$orig.ident, sub@meta.data$Phase))
colnames(stat) <- c("Cluster", "Treatment", "Phase", "Frequency")

g5 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Phase)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  ggtitle(paste(Gene, "Percentage_by_cluster+Treatment", sep="_"))+
  scale_fill_brewer(palette="Dark2")+ 
  facet_grid( ~ Treatment)
g5

#Gene2 in Gene (MKI67+/- cells in IL24+ cells)
Gene2 <- "MKI67"
sample.info <- data.frame(row.names = attr(sub@active.ident, 'names'))
indexGene <- match(Gene2, rownames(sub@assays$RNA@counts))
sample.info$Gene2 = 'Neg'
sample.info$Gene2[which(sub@assays$RNA@counts[rownames(sub@assays$RNA@counts)[indexGene], ] != 0)] <- 'Pos'
sub <- AddMetaData(object = sub, metadata = sample.info)

UMAPPlot(sub, group.by="Gene2",cols = c("#d9d9d9", "#fb8072"))+ggtitle(Gene2) 

stat <- as.data.frame(table(sub@active.ident, sub@meta.data$orig.ident, sub@meta.data$Gene2))
colnames(stat) <- c("Cluster", "Treatment", "Expression", "Frequency")

g6 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Expression)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  ggtitle(paste(Gene2, "Percentage_by_cluster+Treatment", sep="_"))+
  facet_grid( ~ Treatment)+
  scale_fill_brewer(palette="Dark2")
g6

#
#Cell cycle of MKI67+/- cells in IL24+ cells
Idents(sub) <- "Gene2"
sub2 <- subset(sub, idents = "Pos") 
UMAPPlot(sub2,cols = "tomato")+ggtitle("MKI67+ in IL24+")

Idents(sub) <- "Gene2"
sub3 <- subset(sub, idents = "Neg") 
UMAPPlot(sub3,cols = "gray")+ggtitle("MKI67- in IL24+")

Idents(sub2) <- "seurat_clusters"
UMAPPlot(sub2)+ggtitle("MKI67+ in IL24+_Clusters")
UMAPPlot(sub2, split.by="orig.ident")+ggtitle("MKI67+ in IL24+_Clusters")

stat <- as.data.frame(table(sub2@active.ident, sub2@meta.data$orig.ident, sub2@meta.data$Phase))
colnames(stat) <- c("Cluster", "Treatment", "Phase", "Frequency")

g7 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Phase)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  ggtitle(paste(Gene2, "Percentage_by_cluster+Treatment", sep="_"))+
  facet_grid( ~ Treatment)
g7


Idents(sub3) <- "seurat_clusters"
UMAPPlot(sub3)+ggtitle("MKI67- in IL24+_Clusters")
UMAPPlot(sub3, split.by="orig.ident")+ggtitle("MKI67- in IL24+_Clusters")

stat <- as.data.frame(table(sub3@active.ident, sub3@meta.data$orig.ident, sub3@meta.data$Phase))
colnames(stat) <- c("Cluster", "Treatment", "Phase", "Frequency")

g8 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Phase)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  ggtitle(paste(Gene2, "Percentage_by_cluster+Treatment", sep="_"))+ 
  facet_grid( ~ Treatment)
g8

# Analysis_3
#Finding differentially expressed features (cluster biomarkers)#########
#If I need to check original RNA
DefaultAssay(PDX) <- "RNA"

#find all markers of cluster X

cluster0.markers <- FindMarkers(PDX, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 10)

#find markers for every cluster compared to all remaining cells, report only the positive ones
#If I need to change >> @active.assay: chr "integrated" to "RNA"

DefaultAssay(PDX) <- "RNA"

#
PDX.markers <- FindAllMarkers(PDX, min.pct = 0.25, logfc.threshold = 0.25)
top10.PDX.markers <- PDX.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
bottom10.PDX.markers <- PDX.markers %>% group_by(cluster) %>% top_n(-10, avg_logFC)

#GS3

write.csv(PDX.markers, "GS3_intergrate_8_RNA_DEG.csv")
write.csv(top10.PDX.markers, file = "top10_GS3_intergrate_8_DEG.csv")
write.csv(bottom10.PDX.markers, file = "bot10_GS3_intergrate_8_DEG.csv")

#SC31

write.csv(PDX.markers, "SC31_DEG_8_RNA.csv")
write.csv(top10.PDX.markers, file = "top10_SC31_DEG_8_RNA.csv")
write.csv(bottom10.PDX.markers, file = "bot20_SC31_DEG_8_RNA.csv")

#In Cx_DEG_Placebo_vs_E2

#GS3

DefaultAssay(PDX) <- "RNA"
Idents(PDX) <- "seurat_clusters"
UMAPPlot(PDX)
sub <- subset(PDX, idents = "7") 
UMAPPlot(sub,cols = "#d9d9d9")+ggtitle("C0")
Idents(sub) <- "orig.ident"
Marker <- FindAllMarkers(sub, only.pos = T, min.pct = 0.25)
UMAPPlot(sub, cols=c("#0066FF", "#FF99FF"))+ggtitle(paste("C3", "by_treatment", sep="_"))
write.csv(Marker, paste(7, "C7_placebo_vs_E2.csv", sep=""))

#SC31

DefaultAssay(PDX) <- "RNA"
Idents(PDX) <- "seurat_clusters"
UMAPPlot(PDX)
sub <- subset(PDX, idents = "7") 
UMAPPlot(sub,cols = "#d9d9d9")+ggtitle("C0")
Idents(sub) <- "orig.ident"
Marker <- FindAllMarkers(sub, only.pos = T, min.pct = 0.25)
UMAPPlot(sub, cols=c("#FF99FF" , "#0066FF"))+ggtitle(paste("C2", "by_treatment", sep="_"))
write.csv(Marker, paste(7, "C7_placebo_vs_E2.csv", sep=""))

#Heat_Map

#No requested features found in the scale.data slot for the RNA assay.
PDX <- ScaleData(object = PDX, features = rownames(PDX))
#
DoHeatmap(PDX, features = top10.PDX.markers$gene) + NoLegend()
  
+scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
+scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
+scale_fill_gradientn(colors = c("blue", "white", "red"))

# Analysis_4
#UMAP_Pos&Neg plot

Gene <- "IL24"

sample.info <- data.frame(row.names = attr(PDX@active.ident, 'names'))
indexGene <- match(Gene, rownames(PDX@assays$RNA@counts))
sample.info$GENE = 'Neg'
sample.info$GENE[which(PDX@assays$RNA@counts[rownames(PDX@assays$RNA@counts)[indexGene], ] != 0)] <- 'Pos'
PDX <- AddMetaData(object = PDX, metadata = sample.info)
UMAPPlot(PDX, group.by="GENE", cols=c("lightgrey", "tomato"), pt.size = 0.1)+
  ggtitle(Gene) 
UMAPPlot(PDX, group.by="GENE", cols=c("lightgrey", "tomato"), pt.size = 0.1, split.by="orig.ident")+
  ggtitle(Gene) 


#Percentage(GENE_Pos/Neg)_cluster/treatment/cell_cycle_Bar_Graph

Gene <- "IL24"

sample.info <- data.frame(row.names = attr(PDX@active.ident, 'names'))
indexGene <- match(Gene, rownames(PDX@assays$RNA@counts))
sample.info$GENE = 'Neg'
sample.info$GENE[which(PDX@assays$RNA@counts[rownames(PDX@assays$RNA@counts)[indexGene], ] != 0)] <- 'Pos'
PDX <- AddMetaData(object = PDX, metadata = sample.info)

#
Idents(PDX) <- "GENE"
sub <- subset(PDX, idents = "Pos") 
UMAPPlot(sub,cols = "#d9d9d9")+ggtitle("ESR1+")
Idents(sub) <- "orig.ident"

stat <- as.data.frame(table(PDX@active.ident, PDX@meta.data$orig.ident, PDX@meta.data$GENE))
colnames(stat) <- c("Cluster", "Treatment", "Expression", "Frequency")

g4 <- ggplot(data=stat, aes(x=Treatment, y=Frequency, fill=Expression)) +
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  ggtitle(paste(Gene, "Percentage_by_treatment", sep="_"))+
  scale_fill_brewer(palette="Dark2")
g4
write.csv(stat, paste(Gene, ".csv", sep=""))

g5 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Expression)) +
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  ggtitle(paste(Gene, "Percentage_by_cluster+Treatment", sep="_"))+
  scale_fill_brewer(palette="Dark2")+ 
  facet_grid( ~ Treatment)
g5


# Analysis_5
#Percentage2_2genes_Combination

#UMAP

Gene1 <- "IL24"
sample.info <- data.frame(row.names = attr(PDX@active.ident, 'names'))
indexGene <- match(Gene1, rownames(PDX@assays$RNA@counts))
sample.info$GENE1 = 'Neg'
sample.info$GENE1[which(PDX@assays$RNA@counts[rownames(PDX@assays$RNA@counts)[indexGene], ] != 0)] <- 'Pos'
PDX <- AddMetaData(object = PDX, metadata = sample.info)

UMAPPlot(PDX, group.by="GENE1")
UMAPPlot(PDX, group.by="GENE1",cols = c("#d9d9d9", "#fb8072"))+ggtitle(Gene1) 
UMAPPlot(PDX, group.by="GENE1",cols = c("#d9d9d9", "#fb8072"), split.by="orig.ident")+ggtitle(Gene1) 

Gene2 <- "ESR1"
sample.info <- data.frame(row.names = attr(PDX@active.ident, 'names'))
indexGene <- match(Gene2, rownames(PDX@assays$RNA@counts))
sample.info$GENE2 = 'Neg'
sample.info$GENE2[which(PDX@assays$RNA@counts[rownames(PDX@assays$RNA@counts)[indexGene], ] != 0)] <- 'Pos'
PDX <- AddMetaData(object = PDX, metadata = sample.info)

UMAPPlot(PDX, group.by="GENE2")
UMAPPlot(PDX, group.by="GENE2",cols = c("#d9d9d9", "#fb8072"))+ggtitle(Gene2) 
UMAPPlot(PDX, group.by="GENE2",cols = c("#d9d9d9", "#fb8072"), split.by="orig.ident")+ggtitle(Gene1) 

#Bar graph

Gene1 <- "ESR1"
Gene2 <- "IL24"

sample.info <- data.frame(row.names = attr(PDX@active.ident, 'names'))
indexGene <- match(Gene1, rownames(PDX@assays$RNA@counts))
sample.info$Gene1 = 'Neg'
sample.info$Gene1[which(PDX@assays$RNA@counts[rownames(PDX@assays$RNA@counts)[indexGene], ] != 0)] <- 'Pos'
PDX <- AddMetaData(object = PDX, metadata = sample.info)

sample.info <- data.frame(row.names = attr(PDX@active.ident, 'names'))
indexGene <- match(Gene2, rownames(PDX@assays$RNA@counts))
sample.info$Gene2 = 'Neg'
sample.info$Gene2[which(PDX@assays$RNA@counts[rownames(PDX@assays$RNA@counts)[indexGene], ] != 0)] <- 'Pos'
PDX <- AddMetaData(object = PDX, metadata = sample.info)

stat <- as.data.frame(table(PDX$Gene1, PDX$Gene2, PDX$orig.ident))
colnames(stat) <- c("Gene1","Gene2", "Treatment", "Frequency")

ggplot(data=stat, aes(x=Gene1, y=Frequency, fill=Gene2)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  xlab(Gene1)+
  scale_fill_brewer(palette="Dark2")

ggplot(data=stat, aes(x=Gene1, y=Frequency, fill=Gene2)) +
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  labs(x = Gene1, y = "Frequency", fill = Gene2)
xlab(Gene1)+
  scale_fill_brewer(palette="Dark2")

#Percentage2_Bar_Graph_gene_Combination_###separated by treatment##### 

ggplot(data=stat, aes(x=Gene1, y=Frequency, fill=Gene2)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  xlab(Gene1)+
  scale_fill_brewer(palette="Dark2")+
  facet_grid(.~Treatment)

ggplot(data=stat, aes(x=Gene1, y=Frequency, fill=Gene2)) +
  geom_bar(stat="identity", position = "fill")+
  theme_classic()+
  labs(x = Gene1, y = "Frequency", fill = Gene2)+
  xlab(Gene1)+
  ggtitle(paste(Gene2, Gene1, sep="_"))+  
  scale_fill_brewer(palette="Dark2") +
  facet_grid(.~Treatment) 

title <- paste(gene1, gene2, sep="_")
write.csv(stat, paste(title, "_expression.csv", sep=""))  

# Analysis_6
#Cell cycle_by_gene

Gene <- "IL24"

sample.info <- data.frame(row.names = attr(PDX@active.ident, 'names'))
indexGene <- match(Gene, rownames(PDX@assays$RNA@counts))
sample.info$GENE = 'Neg'
sample.info$GENE[which(PDX@assays$RNA@counts[rownames(PDX@assays$RNA@counts)[indexGene], ] != 0)] <- 'Pos'
PDX <- AddMetaData(object = PDX, metadata = sample.info)

stat <- as.data.frame(table(PDX$Phase, PDX$GENE))
colnames(stat) <- c("Phase", "Expression", "Frequency")
ggplot(data=stat, aes(x=Expression, y=Frequency, fill=Phase)) +
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  ggtitle(paste(Gene, "Cell cycle", sep="_"))+
  scale_fill_brewer(palette="Dark2")

write.csv(stat, paste(Gene, "_cell_cycle.csv", sep=""))


Idents(PDX) <- "seurat_clusters"


#Cell cycle+DEG separated by treatment

Gene <- "ESR1"

sample.info <- data.frame(row.names = attr(PDX@active.ident, 'names'))
indexGene <- match(Gene, rownames(PDX@assays$RNA@counts))
sample.info$GENE = 'Neg'
sample.info$GENE[which(PDX@assays$RNA@counts[rownames(PDX@assays$RNA@counts)[indexGene], ] != 0)] <- 'Pos'
PDX <- AddMetaData(object = PDX, metadata = sample.info)

stat <- as.data.frame(table(PDX$Phase, PDX$GENE, PDX$orig.ident))
colnames(stat) <- c("Phase", "Expression", "Treatment", "Frequency")
ggplot(data=stat, aes(x=Expression, y=Frequency, fill=Phase)) +
  geom_bar(stat="identity", position = "fill")+
  theme_classic()+
  ggtitle(paste(Gene, "Cell cycle", sep="_"))+
  scale_fill_brewer(palette="Dark2")+
  facet_grid(.~Treatment)+xlab(Gene)

#cell_number
write.csv(stat, paste(Gene, "_cell_cycle_by_treatment.csv", sep=""))


#Cell_cycle_Bar_Graph

stat <- as.data.frame(table(PDX@active.ident, PDX@meta.data$orig.ident, PDX@meta.data$Phase))
colnames(stat) <- c("Cluster", "Treatment", "Cycle", "Frequency")

g21 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Cycle)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  ggtitle(paste("Frequency_by_cluster", sep="_"))+
  scale_fill_brewer(palette="Dark2")
g21

g22 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Cycle)) +
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  ggtitle(paste("Percentage_by_cluster", sep="_"))+
  scale_fill_brewer(palette="Dark2")
g22

g23 <- ggplot(data=stat, aes(x=Treatment, y=Frequency, fill=Cycle)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  ggtitle(paste("Frequency_by_treatment", sep="_"))+
  scale_fill_brewer(palette="Dark2")
g23

g24 <- ggplot(data=stat, aes(x=Treatment, y=Frequency, fill=Cycle)) +
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  ggtitle(paste("Percentage_by_treatment", sep="_"))+
  scale_fill_brewer(palette="Dark2")
g24

g25 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Cycle)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  ggtitle(paste("Frequency_by_cluster+Treatment", sep="_"))+
  scale_fill_brewer(palette="Dark2")+ 
  facet_grid( ~ Treatment)
g25

g26 <- ggplot(data=stat, aes(x=Cluster, y=Frequency, fill=Cycle)) +
  geom_bar(stat="identity", position = "fill")+
  theme_minimal()+
  ggtitle(paste("Percentage_by_cluster+Treatment", sep="_"))+
  scale_fill_brewer(palette="Dark2")+ 
  facet_grid( ~ Treatment)
g26

CombinePlots(plots = list(g21, g22, g23, g24, g25, g26), ncol=2)

#
write.csv(stat, "CellCycle.csv")

#
stat <- as.data.frame(table(PDX@meta.data$orig.ident, PDX@meta.data$Phase))
write.csv(stat, "CellCycle2.csv")


# Analysis_7
DEG_Analysis

#E2 vs Placebo

Idents(PDX) <- "orig.ident"
UMAPPlot(PDX, cols=c("#FF99FF", "#0066FF"))
DEG_E2vsPlacebo <- FindAllMarkers(PDX, only.pos = T)

Exp <- AverageExpression(PDX, assays = "RNA")
Exp <- Exp[[1]]
Exp <- cbind(row.names(Exp), Exp)
colnames(Exp)[1] <- "gene"

DEG_E2vsPlacebo <- merge(DEG_E2vsPlacebo, Exp)
write.csv(DEG_E2vsPlacebo, "DEG_E2vsPlacebo.csv")

#DEG_Separated_by_GENE

Idents(PDX) <- "GENE"
Marker <- FindAllMarkers(PDX, only.pos = T, min.pct = 0.25)
write.csv(Marker, paste(Gene, "_DEG.csv", sep=""))

#In ESR1+ cell_DEG_Placebo_vs_E2

Idents(PDX) <- "GENE"
sub <- subset(PDX, idents = "Pos") 
UMAPPlot(sub,cols = "#d9d9d9")+ggtitle("ESR1+")
Idents(sub) <- "orig.ident"
Marker <- FindAllMarkers(sub, only.pos = T, min.pct = 0.25)
UMAPPlot(sub, cols=c("#0066FF", "#FF99FF"))+ggtitle(paste(Gene, "+_by_treatment", sep="_"))
write.csv(Marker, paste(Gene, "+_vehicle_vs_E2.csv", sep=""))

#In ESR1- cell_DEG_Placebo_vs_E2

sub <- subset(PDX, idents = "Neg") 
UMAPPlot(sub,cols = "#d9d9d9")+ggtitle("ESR1-")
Idents(sub) <- "orig.ident"
Marker <- FindAllMarkers(sub, only.pos = T, min.pct = 0.25)
UMAPPlot(sub, cols=c("#0066FF", "#FF99FF"))+ggtitle(paste(Gene, "-_by_treatment", sep="_"))
write.csv(Marker, paste(Gene, "-_vehicle_vs_E2.csv", sep=""))

#In Placebo_DEG_ESR1+_vs_ESR1-

Idents(PDX) <- "orig.ident"
sub <- subset(PDX, idents = "1-Placebo") 
UMAPPlot(sub, cols = "#0066FF")+ggtitle("Placebo")
Idents(sub) <- "GENE"
Marker <- FindAllMarkers(sub, only.pos = T, min.pct = 0.25)
UMAPPlot(sub, cols=c("#fb8072","#d9d9d9"))+ggtitle("Placebo_ESR1+/-")
write.csv(Marker, paste(Gene, "_DEG_in_Placebo.csv", sep=""))

#In E2_DEG_ESR1+_vs_ESR1-

sub <- subset(PDX, idents = "2-E2") 
UMAPPlot(sub, cols="#FF99FF")+ggtitle("E2")
Idents(sub) <- "GENE"
Marker <- FindAllMarkers(sub, only.pos = T, min.pct = 0.25)
UMAPPlot(sub, cols=c("#d9d9d9", "#fb8072"))+ggtitle("E2_ESR1+/-")
write.csv(Marker, paste(Gene, "_DEG_in_E2.csv", sep=""))


remove(sub)
gc()
gc()

Idents(PDX) <- "seurat_clusters"
