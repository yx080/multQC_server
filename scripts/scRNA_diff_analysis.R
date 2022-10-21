args <- commandArgs(trailingOnly = TRUE)
rds=args[1]
diff_results=args[2]
diff_path=args[3]
p.feature=args[4]
rna_ump_markergene=args[5]
rna_heatmap=args[6]
bc_mtx=args[7]
i.matrix=args[8]

#rds="analysis/wt/rna_QC/wt.scRNA.qc.rds"
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
#diff_results="analysis/wt/rna_QC/diff_analysis/wt.diff_results.csv"
#diff_path="analysis/wt/rna_QC/diff_analysis"
#p.feature="Gapdh,Plcb1"
#rna_ump_markergene="analysis/f11/rna_QC/wt.rna_ump_markergene.png"
#rna_heatmap="analysis/wt/rna_QC/wt.rna_heatmap.png"
#bc_mtx="analysis/f11/joint_cell_calling/wt.barcodes.csv"
#i.matrix="data/matrix/wt.h5"

#p.feature<-as.character(toupper(unlist(strsplit(p.feature,split = ","))))
p.feature<-as.character(unlist(strsplit(p.feature,split = ",")))

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

#matx_new2<-read.csv(bc_mtx)
#counts <- Read10X_h5(i.matrix)
subset_rna<-readRDS(rds)
DefaultAssay(subset_rna) <- "RNA"

all.markers <- FindAllMarkers(object = subset_rna)
head(x = all.markers)
dir.create(paste(diff_path,sep=""))
write.csv(all.markers,diff_results,quote = F)

cluster_rna_num<-length(levels(subset_rna))-1

all.markers$gene<-row.names(all.markers)
all.markers$gene<-gsub("\\..*", "", all.markers$gene)

#dir.create("multiome_QC_results/rna_QC/diff_analysis")
rna_de.analysis<-function(i){
  plot1 <- VlnPlot(subset_rna, features = c(filter(all.markers,cluster==i)[1,7]))+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6),axis.title.y = element_text(color = "grey20", size = 4)) 
  
  plot2 <- FeaturePlot(
    object = subset_rna,
    features = filter(all.markers,cluster==i)[1,7],
    pt.size = 0.1,
    max.cutoff = 'q95'
  )
  png(file=paste(diff_path,"/rna_diff_analysis_", i,".png", sep=""), width = 6, height = 3, unit="in", res=200)
  print(plot1+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6),axis.title.y = element_text(color = "grey20", size = 4)) | plot2+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6),axis.title.y = element_text(color = "grey20", size = 4)))
  dev.off()
}

#for (i in 0:cluster_rna_num) {
#  rna_de.analysis(i)
#}

png(file=rna_ump_markergene, width = 10, height = 5, unit="in", res=200)
FeaturePlot(
  object = subset_rna,
  features = p.feature,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
dev.off()


heatmp_feature <-""
for (i in 0:cluster_rna_num) {
  heatmp_feature <- append(heatmp_feature,c(dplyr::filter(all.markers,cluster==i)[c(1:5),7]))
}
heatmp_feature<-heatmp_feature[-1]

subset_rna$groups <- sample(c("group1", "group2"), size = ncol(subset_rna), replace = TRUE)
features <- heatmp_feature

#RidgePlot(subset_rna, features = features, ncol = 2)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6),axis.title.y = element_text(color = "grey20", size = 4))
png(file=rna_heatmap, width = 6, height = 6, unit="in", res=200)
DoHeatmap(subset(subset_rna, downsample = 100), features = features, size = 2) 
dev.off()

