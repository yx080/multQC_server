args <- commandArgs( trailingOnly = TRUE )
rds=args[1]
diff_path=args[2]
atac_umap_marker_gene=args[3]
p.feature=args[4]
#motif_path=args[5]
p.feature<-as.character(unlist(strsplit(p.feature,split = ",")))
genome=args[5]

#rds="analysis/wt/atac_QC/wt.scATAC.qc.rds"
#matx_file="analysis/f9/joint_cell_calling/f9.barcodes.csv"
#peak_target="analysis/f9/atac_QC/f9.peak_target.png"
#tss_insertsize="analysis/f9/atac_QC/f9.tss_insertsize.png"
#atac_qc_violin="analysis/f9/atac_QC/f9.atac_qc_violin.png"
#atac_depthCor="analysis/f9/atac_QC/f9.atac_depthCor.png"
#atac_umap_cluster="analysis/f9/atac_QC/f9.atac_umap_cluster.png"
#atac_coordinates="analysis/f9/atac_QC/f9.atac_coordinates.csv" 
#atac_umap_byvennbarcode="analysis/f9/atac_QC/f9.atac_umap_byvennbarcode.png"
#atac_umap_depth="analysis/f9/atac_QC/f9.atac_umap_depth.png"
#diff_path="analysis/wt/atac_QC/diff_analysis"
#atac_umap_marker_gene="analysis/wt/atac_QC/wt.atac_umap_marker_gene.png" 
#p.feature="Top2a,Cenpe"
#motif_path="analysis/f9/atac_QC/motif_analysis"
#genome='mm10'

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
#library(future)


libs<-""
ifelse(genome=="mm10", {libs[1] <- "EnsDb.Mmusculus.v79" ; libs[2]<-"BSgenome.Mmusculus.UCSC.mm10";},
       {libs[1]<-"EnsDb.Hsapiens.v86" ; libs[2]<-"BSgenome.Hsapiens.UCSC.hg38";})

for (i in libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}

subset_atac <- readRDS(rds)
DefaultAssay(subset_atac) <- "ATAC"
#matx_new <- read.csv(matx_file)


#plan("multiprocess", workers = 8)
cluster_num<-length(levels(subset_atac))-1

#dir.create(paste(sampleName,"/atac_QC/diff_analysis",sep = ""))
dir.create(diff_path)
de_analysis<-function(i) {
  da_peaks <- FindMarkers(
    object = subset_atac,
    ident.1 = i, 
    ident.2 = NULL,
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nCount_ATAC'
  )
  open_c0 <- rownames(da_peaks[da_peaks$avg_log2FC > 0, ])
  if(nrow(da_peaks)==0){
    text_print<-"no_differential_peaks_found"
    write.csv(text_print,file = paste(diff_path, "/atac_diff_analysis_",i,".csv",sep = ""),quote = F)
  }else{
    closest_genes_c0<- ClosestFeature(subset_atac, regions = open_c0)
    write.csv(closest_genes_c0[,c(2,4,5,6,7)],file = paste(diff_path, "/atac_diff_analysis_",i,".csv",sep = ""),quote = F)
    
    #setdiff(0:cluster_num,0)
    
    plot1 <- VlnPlot(
      object = subset_atac,
      features = rownames(da_peaks)[1],
      pt.size = 0.1,
      idents = c(0:cluster_num)
    )
    plot2 <- FeaturePlot(
      object = subset_atac,
      features = rownames(da_peaks)[1],
      pt.size = 0.1,
      max.cutoff = 'q95'
    )
    png(file=paste(diff_path,"/atac_diff_analysis_", i,".png", sep=""), width = 4, height = 2, unit="in", res=200)
    print(plot1+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6),axis.title.y = element_text(color = "grey20", size = 4)) | plot2+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6),axis.title.y = element_text(color = "grey20", size = 4)))
    dev.off()
  }
}

for (i in 0:cluster_num) {
  de_analysis(i)
}

gene.activities <- GeneActivity(subset_atac)

# add the gene activity matrix to the Seurat object as a new assay
subset_atac[['nRNA']] <- CreateAssayObject(counts = gene.activities)
subset_atac <- NormalizeData(
  object = subset_atac,
  assay = 'nRNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(subset_atac$nCount_RNA)
)

DefaultAssay(subset_atac) <- 'nRNA'

#parameter
#test.feature<-c('Psca','Krt5',"Ar","Cd14","Krt8","Cd44")

#p.feature<-as.character(toupper(unlist(strsplit(p.feature,split = ","))))
p.feature<-as.character(unlist(strsplit(p.feature,split = ",")))

rownames(subset_atac)

png(file=atac_umap_marker_gene, width = 9, height = 5.8, unit="in", res=200)
FeaturePlot(
  object = subset_atac,
  features = p.feature,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3)
dev.off()

# DefaultAssay(subset_atac) <- "ATAC"

# #parameter
# #p.region<-c("Psca", "Krt5")

# #png(file=paste("multiome_QC_results/atac_QC/","atac_genomic_region.png", sep=""), width = 6, height = 5, unit="in", res=200)
# #CoveragePlot(
# #  object = subset_atac,
# #  region = p.region,
# #  extend.upstream = 1000,
# #  extend.downstream = 1000,
# #  ncol = 1
# #) + theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 6),axis.text.x = element_text(color = "grey20", size = 6),axis.title.x = element_text(color = "grey20", size = 6),title =element_text(size=6),axis.title.y = element_text(color = "grey20", size = 6))
# #dev.off()

# # extract position frequency matrices for the motifs
# pwm <- getMatrixSet(
#   x = JASPAR2020,
#   opts = list(species = 9606, all_versions = FALSE)
# )

# # add motif information
# # DimPlot(subset_atac, label = TRUE)
# if(genome == "mm10"){
#   main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
# } else {  
#   main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
# }

# keep.peaks <- which(as.character(seqnames(granges(subset_atac))) %in% main.chroms)
# subset_atac[["ATAC"]] <- subset(subset_atac[["ATAC"]], features = rownames(subset_atac[["ATAC"]])[keep.peaks])

# if(genome == "mm10"){
#   subset_atac <- AddMotifs(subset_atac, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pwm)
# } else {  
#   subset_atac <- AddMotifs(subset_atac, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)
# }


# if(genome == "mm10"){
#   subset_atac <- RunChromVAR(
#     object = subset_atac,
#     genome = BSgenome.Mmusculus.UCSC.mm10
#   )
# } else {  
#   subset_atac <- RunChromVAR(
#     object = subset_atac,
#     genome = BSgenome.Hsapiens.UCSC.hg38
#   )
# }



# DefaultAssay(subset_atac) <- 'chromvar'
# #dir.create(paste(sampleName,"/atac_QC/motif_analysis",sep = ""))
# dir.create(motif_path)
# motif_analysis<-function(k){
#   differential.activity <- FindMarkers(
#     object = subset_atac,
#     ident.1 = k,
#     ident.2 = NULL,
#     only.pos = TRUE,
#     mean.fxn = rowMeans,
#     fc.name = "avg_diff"
#   )
  
#   write.csv(differential.activity,paste(motif_path,"/motif_analysis_cluster",k,".csv", sep=""),quote = F)
  
#   png(file=paste(motif_path,"/motif_analysis_cluster",k,".png", sep=""), width = 8, height = 4, unit="in", res=200)    
#   print(MotifPlot(
#     object = subset_atac,
#     motifs = head(rownames(differential.activity)),
#     assay = 'ATAC'
#   ))
#   dev.off()
# }

# for (i in 0:cluster_num) {
#   motif_analysis(i)
# }
