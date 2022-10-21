args <- commandArgs( trailingOnly = TRUE )
rds=args[1]
matx_file=args[2]
peak_target=args[3]
tss_insertsize=args[4]
atac_qc_violin=args[5]
atac_depthCor=args[6]
atac_umap_cluster=args[7]
atac_coordinates=args[8]
atac_umap_byvennbarcode=args[9]
atac_umap_depth=args[10]
#diff_path=args[11]
#atac_umap_marker_gene=args[12]
#p.feature=args[13]
#motif_path=args[14]
#p.feature<-as.character(unlist(strsplit(p.feature,split = ",")))
genome=args[11]
atac_outrds=args[12]
#rds="analysis/f9/atac_QC/f9.scATAC.rds"
#matx_file="analysis/f9/joint_cell_calling/f9.barcodes.csv"
#peak_target="analysis/f9/atac_QC/f9.peak_target.png"
#tss_insertsize="analysis/f9/atac_QC/f9.tss_insertsize.png"
#atac_qc_violin="analysis/f9/atac_QC/f9.atac_qc_violin.png"
#atac_depthCor="analysis/f9/atac_QC/f9.atac_depthCor.png"
#atac_umap_cluster="analysis/f9/atac_QC/f9.atac_umap_cluster.png"
#atac_coordinates="analysis/f9/atac_QC/f9.atac_coordinates.csv" 
#atac_umap_byvennbarcode="analysis/f9/atac_QC/f9.atac_umap_byvennbarcode.png"
#atac_umap_depth="analysis/f9/atac_QC/f9.atac_umap_depth.png"
#diff_path="analysis/f9/atac_QC/diff_analysis"
#atac_umap_marker_gene="analysis/f9/atac_QC/f9.atac_umap_marker_gene.png" 
#p.feature="Top2a,Cenpe"
#motif_path="analysis/f9/atac_QC/motif_analysis"
#genome='hg38'

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

libs<-""
ifelse(genome=="mm10", {libs[1] <- "EnsDb.Mmusculus.v79" ; libs[2]<-"BSgenome.Mmusculus.UCSC.mm10";},
       {libs[1]<-"EnsDb.Hsapiens.v86" ; libs[2]<-"BSgenome.Hsapiens.UCSC.hg38";})

for (i in libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}

subset_atac <- readRDS(file = rds)
DefaultAssay(subset_atac) <- "ATAC"
matx_new <- read.csv(matx_file)

peaktarget_10x<-ggplot(matx_new, aes(x=atac_fragments, y=frip,colour=iscell)) +
  geom_point(size=0.01) 
#p3 + scale_x_log10()
matx_new$ataccells<-ifelse(matx_new$cells=='atac_cells',"yes","no")
matx_new$ataccells<-matx_new$ataccells %>% replace_na('no')
peaktarget_custom<-ggplot(matx_new, aes(x=atac_fragments, y=frip,colour=ataccells)) +
  geom_point(size=0.01) 

#dir.create(paste(sampleName,"/atac_QC",sep=""))
png(file=peak_target, width = 10, height = 5, unit="in", res=200)
peaktarget_10x+scale_x_log10()+ theme_minimal(base_size = 10)+theme(legend.position = "none") | peaktarget_custom +scale_x_log10()+ theme_minimal(base_size = 10)
dev.off()

subset_atac <- NucleosomeSignal(object = subset_atac)
subset_atac<- TSSEnrichment(object = subset_atac, fast = FALSE)
tss<-TSSPlot(subset_atac) + NoLegend()+ theme_minimal(base_size = 10) + theme(legend.position = "none")+labs(title="")
subset_atac$nucleosome_signal
insertsize<-FragmentHistogram(object = subset_atac,region = 'chr1-1-10000000')+ theme_minimal(base_size = 10) + theme(legend.position = "none")+labs(title="")
png(file=tss_insertsize, width = 10, height = 5, unit="in", res=200)
tss | insertsize
dev.off()

v1<-VlnPlot(
  object = subset_atac,
  features = c('nCount_ATAC'),
  pt.size = 0,ncol = 3)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 8),axis.text.x = element_text(color = "white", size = 8),axis.title.x = element_text(color = "white", size = 8),legend.title=element_text(size=6),title =element_text(size=8))
v2<-VlnPlot(
  object = subset_atac,
  features = c('nFeature_ATAC'),
  pt.size = 0,ncol = 3)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 8),axis.text.x = element_text(color = "white", size = 8),axis.title.x = element_text(color = "white", size = 8),legend.title=element_text(size=6),title =element_text(size=8))
v3<-VlnPlot(
  object = subset_atac,
  features = c('TSS.enrichment'),
  pt.size = 0,ncol = 3)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 8),axis.text.x = element_text(color = "white", size = 8),axis.title.x = element_text(color = "white", size = 8),legend.title=element_text(size=6),title =element_text(size=8))
v4<-VlnPlot(
  object = subset_atac,
  features = c('nucleosome_signal'),
  pt.size = 0,ncol = 3)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 8),axis.text.x = element_text(color = "white", size = 8),axis.title.x = element_text(color = "white", size = 8),title =element_text(size=8))

png(file=atac_qc_violin, width = 8, height = 4, unit="in", res=200)
grid.arrange(v1,v2,v3,v4,ncol = 4)
dev.off()

subset_atac <- RunTFIDF(subset_atac)
subset_atac <- FindTopFeatures(subset_atac, min.cutoff = 'q0')
subset_atac <- RunSVD(subset_atac)

png(file=atac_depthCor, width = 6, height = 5, unit="in", res=200)
DepthCor(subset_atac)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 7),axis.text.x = element_text(color = "grey20", size = 7),axis.title.x = element_text(color = "grey20", size = 7),title =element_text(size=7),axis.title.y = element_text(color = "grey20", size = 7))
dev.off()

subset_atac <- RunUMAP(object = subset_atac, reduction = 'lsi', dims = 2:30)
subset_atac <- FindNeighbors(object = subset_atac, reduction = 'lsi', dims = 2:30)
subset_atac <- FindClusters(object = subset_atac, verbose = FALSE, algorithm = 3)

png(file=atac_umap_cluster, width = 6, height = 5, unit="in", res=200)
DimPlot(object = subset_atac, label = TRUE,pt.size =0.001)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 7),axis.text.x = element_text(color = "grey20", size = 7),axis.title.x = element_text(color = "grey20", size = 7),title =element_text(size=7),axis.title.y = element_text(color = "grey20", size = 7))
dev.off()

coordinates.atac<-as.data.frame(subset_atac[["umap"]]@cell.embeddings)
write.csv(coordinates.atac,file=atac_coordinates,quote = F)
coor2<-cbind(coordinates.atac,subset_atac$nCount_RNA,subset_atac$nCount_ATAC)

cellsinboth<-dplyr::filter(matx_new,barcodes_annotation=='cells_in_both')
add_anno<-as.data.frame(cbind(cellsinboth$barcode,cellsinboth$barcodes_annotation))
colnames(add_anno)<-c("barcode","annotation")
coor2$barcode<-row.names(coor2)
coor2_new<-left_join(coor2,add_anno)
coor2_new[is.na(coor2_new)] <- "atac_cells"

u1<-ggplot(coor2_new,aes(x=UMAP_1,y=UMAP_2,color=annotation))+ geom_point(size=0.01) + scale_color_manual(values = c("cells_in_both" = "orange","atac_cells"="blue")) +theme_classic()+theme(legend.position = "right", axis.text.y = element_text(color = "grey20", size = 7),legend.key.size = unit(0.2, 'cm'),axis.text.x = element_text(color = "grey20", size = 7),legend.title = element_text(size=7),legend.text = element_text(size=7),axis.title.x = element_text(color = "grey20", size = 7),title =element_text(size=7),axis.title.y = element_text(color = "grey20", size = 7))
png(file=atac_umap_byvennbarcode, width = 6, height = 5, unit="in", res=200)
u1
dev.off()

coor2$Depth<-log10(subset_atac$nCount_ATAC)

png(file=atac_umap_depth, width = 6, height = 5, unit="in", res=200)
ggplot(coor2,aes(x=UMAP_1,y=UMAP_2,color=Depth))+theme_classic()+ geom_point(size=0.01)+ scale_color_gradientn(colours = viridis(5))+ theme(legend.position = "right", axis.text.y = element_text(color = "grey20", size = 7),legend.key.size = unit(0.2, 'cm'),axis.text.x = element_text(color = "grey20", size = 7),legend.title = element_text(size=7),legend.text = element_text(size=7),axis.title.x = element_text(color = "grey20", size = 7),title =element_text(size=7),axis.title.y = element_text(color = "grey20", size = 7))
dev.off()

saveRDS(subset_atac,atac_outrds)