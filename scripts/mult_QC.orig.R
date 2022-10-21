args <- commandArgs(trailingOnly = TRUE)
rds=args[1]
jointumap=args[2]
mult_coordinates=args[3]
jointumap_rna=args[4]
jointumap_atac=args[5]
genometrack_path=args[6]
stats_multtb=args[7]
genome=args[8]
p.feature=args[9]
bc_mtx=args[10]
mult_newrds=args[11]
#rds="analysis/f11/mult_QC/f11.mult.rds"
#jointumap="analysis/f11/mult_QC/f11.jointumap.png" 
#mult_coordinates="analysis/f11/mult_QC/f11.mult_coordinates.csv"
#jointumap_rna="analysis/f11/mult_QC/f11.jointumap_rna.png" 
#jointumap_atac="analysis/f11/mult_QC/f11.jointumap_atac.png" 
#genometrack_path="analysis/f11/mult_QC" 
#stats_multtb="analysis/f11/mult_QC/f11.stats_mult.csv" 
#genome="hg38"
#p.feature="Top2a,Cenpe"
#bc_mtx="analysis/f11/joint_cell_calling/f11.barcodes.csv"

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

p.feature1<-as.character(unlist(strsplit(p.feature,split = ",")))
p.feature<-as.character(toupper(unlist(strsplit(p.feature,split = ","))))

libs<-""
ifelse(genome=="mm10", {libs[1] <- "EnsDb.Mmusculus.v79" ; libs[2]<-"BSgenome.Mmusculus.UCSC.mm10";},
       {libs[1]<-"EnsDb.Hsapiens.v86" ; libs[2]<-"BSgenome.Hsapiens.UCSC.hg38";})

for (i in libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}

subset_mult<-readRDS(rds)
subset_mult

DefaultAssay(subset_mult) <- "RNA"
subset_mult <- SCTransform(subset_mult)
subset_mult <- RunPCA(subset_mult)

DefaultAssay(subset_mult) <- "ATAC"
subset_mult <- FindTopFeatures(subset_mult, min.cutoff = 5)
subset_mult <- RunTFIDF(subset_mult)
subset_mult <- RunSVD(subset_mult)

subset_mult <- RunUMAP(
  object = subset_mult,
  assay = "ATAC",
  dims = 2:30
)

subset_mult <- FindNeighbors(object = subset_mult, reduction = 'lsi', dims = 2:30)
subset_mult <- FindClusters(object = subset_mult, verbose = FALSE, algorithm = 3)

subset_mult <- FindMultiModalNeighbors(
  object = subset_mult,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  weighted.nn.name = "weighted.nn",
  verbose = TRUE
)

subset_mult <- RunUMAP(
  object = subset_mult,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

#dir.create(paste(sampleName,"/mult",sep = ""))
png(file=jointumap, width = 5, height = 4, unit="in", res=200)
DimPlot(subset_mult, label = TRUE, repel = TRUE, reduction = "umap",pt.size = 0.01) + NoLegend()+ theme(legend.position = "right",legend.text = element_text(color = "grey20", size = 7), axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 7),axis.title.x = element_text(color = "grey20", size = 7),title =element_text(size=6),axis.title.y = element_text(color = "grey20", size = 7))
dev.off()

coordinates.mult<-as.data.frame(subset_mult[["umap"]]@cell.embeddings)
write.csv(coordinates.mult,mult_coordinates,quote = F)

DefaultAssay(subset_mult) <- "RNA"

png(file=jointumap_rna, width = 10, height = 5, unit="in", res=200)
FeaturePlot(
  object = subset_mult,
  features = p.feature1,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3)
dev.off()

DefaultAssay(subset_mult) <- "ATAC"
gene.activities <- GeneActivity(subset_mult)

# add the gene activity matrix to the Seurat object as a new assay
subset_mult[['activity']] <- CreateAssayObject(counts = gene.activities)
subset_mult <- NormalizeData(
  object = subset_mult,
  assay = 'activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(subset_mult$nCount_RNA)
)


DefaultAssay(subset_mult) <- 'activity'
png(file=jointumap_atac, width = 10, height = 5, unit="in", res=200)
FeaturePlot(
  object = subset_mult,
  features = p.feature1,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
dev.off()

saveRDS(subset_mult, file=mult_newrds)

