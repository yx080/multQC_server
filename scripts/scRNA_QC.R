args <- commandArgs(trailingOnly = TRUE)
rds=args[1]
rna_qc_violin=args[2]
feature_cor=args[3]
rna_highly_variable_gene=args[4]
rna_umap_cluster=args[5]
rna_coordinates=args[6]
rna_umap_byvennbarcode=args[7]
rna_umap_depth=args[8]
rna_umap_mt_percent=args[9]
#rna_ump_markergene=args[10]
bc_mtx=args[10]
i.matrix=args[11]
rna_outrds=args[12]
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
#i.matrix="data/matrix/f11.h5"

#feature<-as.character(toupper(unlist(strsplit(p.feature,split = ","))))

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
#library(motifmatchr)
#library(JASPAR2020)
#library(TFBSTools)
set.seed(1234)

matx_new2<-read.csv(bc_mtx)
counts <- Read10X_h5(i.matrix)
subset_rna<-readRDS(rds)
DefaultAssay(subset_rna) <- "RNA"

subset_rna[["percent.mt"]] <- PercentageFeatureSet(subset_rna, pattern = "^MT-")
r1<-VlnPlot(subset_rna, features = c("nFeature_RNA"),ncol=3,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 6),axis.text.x = element_text(color = "white", size = 6),axis.title.x = element_text(color = "white", size = 6),title =element_text(size=6))
r2<-VlnPlot(subset_rna, features = c("nCount_RNA"),ncol=3,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 6),axis.text.x = element_text(color = "white", size = 6),axis.title.x = element_text(color = "white", size = 6),title =element_text(size=6))
r3<-VlnPlot(subset_rna, features = c("percent.mt"),ncol=3,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 6),axis.text.x = element_text(color = "white", size = 6),axis.title.x = element_text(color = "white", size = 6),title =element_text(size=6))
png(file=rna_qc_violin, width = 5, height = 3, unit="in", res=200) 
grid.arrange(r1,r2,r3,ncol = 4)
dev.off()

plot1 <- FeatureScatter(subset_rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(subset_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(file=feature_cor, width = 5, height = 3, unit="in", res=200) 
plot1 + theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 6),axis.text.x = element_text(color = "grey20", size = 6),axis.title.x = element_text(color = "grey20", size = 6),title =element_text(size=6)) | plot2+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 6),axis.text.x = element_text(color = "grey20", size = 6),axis.title.x = element_text(color = "grey20", size = 6),title =element_text(size=6))
dev.off()

subset_rna <- FindVariableFeatures(subset_rna, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(subset_rna), 10)
plot1 <- VariableFeaturePlot(subset_rna,pt.size = 0.1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(file=rna_highly_variable_gene, width = 4.5, height = 3, unit="in", res=200) 
plot2+ theme(legend.position = "top", text=element_text(size=7),axis.text.y = element_text(color = "grey20", size = 7),axis.text.x = element_text(color = "grey20", size = 7),axis.title.x = element_text(color = "grey20", size = 7),title =element_text(size=7))
dev.off()


subset_rna <- NormalizeData(subset_rna, normalization.method = "LogNormalize", scale.factor = 10000)
subset_rna <- ScaleData(subset_rna)
#subset_rna <- RunPCA(subset_rna, features = VariableFeatures(object = subset_rna))
subset_rna<- FindVariableFeatures(object = subset_rna)
subset_rna<- RunPCA(subset_rna, features = VariableFeatures(object = subset_rna) )

subset_rna <- FindNeighbors(subset_rna, dims = 1:6)
subset_rna <- FindClusters(subset_rna, resolution = 0.5)
subset_rna <- RunUMAP(subset_rna, dims = 1:6)
png(file=rna_umap_cluster, width = 6, height = 5, unit="in", res=200)
DimPlot(subset_rna, reduction = "umap",label = TRUE,pt.size =0.001)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 7),axis.text.x = element_text(color = "grey20", size = 7),axis.title.x = element_text(color = "grey20", size = 7),title =element_text(size=7),axis.title.y = element_text(color = "grey20", size = 7))
dev.off()

coordinates.rna<-as.data.frame(subset_rna[["umap"]]@cell.embeddings)
write.csv(coordinates.rna,rna_coordinates,quote = F)
coor3<-cbind(coordinates.rna, subset_rna$nCount_RNA)
#coor3<-cbind(coordinates.rna, subset_rna$nCount_RNA, subset_rna$nCount_ATAC)
coor3$barcode<-row.names(coor3)
cellsinboth<-dplyr::filter(matx_new2,barcodes_annotation=='cells_in_both')
add_anno<-as.data.frame(cbind(cellsinboth$barcode,cellsinboth$barcodes_annotation))
colnames(add_anno)<-c("barcode","annotation")
coor3_new<-left_join(coor3,add_anno)
coor3_new[is.na(coor3_new)] <- "rna_cells"
#coor2_new$annotation<-as.factor(coor2_new$annotation)

u2<-ggplot(coor3_new,aes(x=UMAP_1,y=UMAP_2,color=annotation))+ geom_point(size=0.01) + scale_color_manual(values = c("cells_in_both" = "orange","rna_cells"="blue")) +theme_classic()+theme(legend.position = "right", axis.text.y = element_text(color = "grey20", size = 7),legend.key.size = unit(0.2, 'cm'),axis.text.x = element_text(color = "grey20", size = 7),legend.title = element_text(size=7),legend.text = element_text(size=7),axis.title.x = element_text(color = "grey20", size = 7),title =element_text(size=7),axis.title.y = element_text(color = "grey20", size = 7))
png(file=rna_umap_byvennbarcode, width = 6, height = 5, unit="in", res=200)
u2
dev.off()

coor3$Depth<-log10(subset_rna$nCount_RNA)
u.rna.1<-ggplot(coor3,aes(x=UMAP_1,y=UMAP_2,color=Depth))+theme_classic()+ geom_point(size=0.01)+ scale_color_gradientn(colours = viridis(5))+ theme(legend.position = "right", axis.text.y = element_text(color = "grey20", size = 7),legend.key.size = unit(0.2, 'cm'),axis.text.x = element_text(color = "grey20", size = 7),legend.title = element_text(size=7),legend.text = element_text(size=7),axis.title.x = element_text(color = "grey20", size = 7),title =element_text(size=7),axis.title.y = element_text(color = "grey20", size = 7))
png(file=rna_umap_depth, width = 6, height = 5, unit="in", res=200)
u.rna.1
dev.off()

mt.rna<-as.data.frame(subset_rna$percent.mt)
mt.rna$barcode<-row.names(mt.rna)
coor3<-left_join(coor3,mt.rna)
coor3$mt.reads.percentage<-coor3$`subset_rna$percent.mt`
u.rna.2<-ggplot(coor3,aes(x=UMAP_1,y=UMAP_2,color=log(mt.reads.percentage)))+theme_classic()+ geom_point(size=0.01)+ scale_color_gradientn(colours = viridis(5))+ theme(legend.position = "right", axis.text.y = element_text(color = "grey20", size = 7),legend.key.size = unit(0.2, 'cm'),axis.text.x = element_text(color = "grey20", size = 4),legend.title = element_text(size=7),legend.text = element_text(size=7),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=7),axis.title.y = element_text(color = "grey20", size = 7))
png(file=rna_umap_mt_percent, width = 6, height = 5, unit="in", res=200)
u.rna.2
dev.off()

saveRDS(subset_rna,rna_outrds)