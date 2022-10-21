args <- commandArgs(trailingOnly = TRUE)
rds=args[1]
rna_umap_cell_cycle=args[2]
rna_cell_cycle=args[3]
bc_mtx=args[4]
i.matrix=args[5]

#rds="analysis/f11/rna_QC/f11.scRNA.rds"
#rna_qc_violin="analysis/f11/rna_QC/f11.rna_qc_violin.png"
#feature_cor="analysis/f11/rna_QC/f11.feature_cor.png" 
#rna_highly_variable_gene="analysis/f11/rna_QC/f11.rna_highly_variable_gene.png"
#rna_umap_cluster="analysis/f11/rna_QC/f11.rna_umap_cluster.png" 
#rna_coordinates="analysis/f11/rna_QC/f11.rna_coordinates.csv"
#rna_umap_byvennbarcode="analysis/f11/rna_QC/f11.rna_umap_byvennbarcode.png" 
#rna_umap_depth="analysis/f11/rna_QC/f11.rna_umap_depth.png" 
#rna_umap_mt_percent="analysis/f11/rna_QC/f11.rna_umap_mt_percent.png"
#rna_umap_cell_cycle="analysis/f11/rna_QC/f11.rna_umap_cell_cycle.png" 
#rna_cell_cycle="analysis/f11/rna_QC/f11.rna_cell_cycle.png" 
#diff_results="analysis/f11/rna_QC/diff_analysis/f11.diff_results.csv"
#diff_path="analysis/f11/rna_QC/diff_analysis"
#p.feature="Top2a,Nlgn1,Cenpe"
#rna_ump_markergene="analysis/f11/rna_QC/f11.rna_ump_markergene.png"
#rna_heatmap="analysis/f11/rna_QC/f11.rna_heatmap.png"
#bc_mtx="analysis/f11/joint_cell_calling/f11.barcodes.csv"
#i.matrix="data/matrix/wt.h5"

#p.feature<-as.character(toupper(unlist(strsplit(p.feature,split = ","))))
#p.feature<-as.character(unlist(strsplit(p.feature,split = ",")))

library(ggplot2)
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(patchwork)
library(dplyr)
library(ggVennDiagram)
library(tidyr)
library(gridExtra)
library(viridis)
library(scales)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
set.seed(1234)

matx_new2<-read.csv(bc_mtx)
counts <- Read10X_h5(i.matrix)
subset_rna<-readRDS(rds)
DefaultAssay(subset_rna) <- "RNA"


ct_rna<-counts[["Gene Expression"]]
ct_rna<-GetAssayData(object = subset_rna, slot = "counts")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

subset_rna <- CellCycleScoring(subset_rna, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
info_rna<-(subset_rna[[]])

png(file=rna_umap_cell_cycle, width = 6, height = 5, unit="in", res=200)
DimPlot(subset_rna,pt.size = 0.01)+ theme(legend.position = "right",legend.text = element_text(color = "grey20", size = 7), axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 7),axis.title.x = element_text(color = "grey20", size = 7),title =element_text(size=6),axis.title.y = element_text(color = "grey20", size = 7))
dev.off()

subset_rna <- FindNeighbors(subset_rna, dims = 1:6)
subset_rna <- FindClusters(subset_rna, resolution = 0.5)
subset_rna <- RunUMAP(subset_rna, dims = 1:6)
rna.cluster<-as.data.frame(subset_rna@active.ident)
rna.cluster$barcode<-rownames(rna.cluster)
row.names(info_rna)->info_rna$barcode
rna.phase<-as.data.frame(cbind(info_rna$barcode,info_rna$Phase))
colnames(rna.phase)<-c("barcode","phase")
rna.id.phase<-merge(rna.cluster,rna.phase)
rna.id.phase$clusters<-rna.id.phase$`subset_rna@active.ident`
png(file=rna_cell_cycle, width = 6, height = 4, unit="in", res=200)
ggplot(rna.id.phase, aes(x = clusters)) + geom_bar(aes(fill = phase), position = 'fill') +
  scale_y_continuous(labels = percent_format()) + theme_bw()+theme_minimal(base_size = 7)
dev.off()

