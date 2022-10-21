args <- commandArgs(trailingOnly = TRUE)
rds=args[1]
#jointumap=args[2]
#mult_coordinates=args[3]
#jointumap_rna=args[4]
#jointumap_atac=args[5]
genometrack_path=args[2]
stats_multtb=args[3]
genome=args[4]
p.feature=args[5]
bc_mtx=args[7]
mult_newrds=args[6]
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

if(genome == "mm10"){
  main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
} else {  
  main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
}


keep.peaks <- which(as.character(seqnames(granges(subset_mult[["ATAC"]]))) == main.chroms)
subset_mult[["ATAC"]] <- subset(subset_mult[["ATAC"]], features = rownames(subset_mult[["ATAC"]])[keep.peaks])

DefaultAssay(subset_mult) <- "ATAC"

if(genome == "mm10"){
  subset_mult <- RegionStats(subset_mult, genome = BSgenome.Mmusculus.UCSC.mm10)
} else {  
  subset_mult <- RegionStats(subset_mult, genome = BSgenome.Hsapiens.UCSC.hg38)
}


subset_mult <- LinkPeaks(
  object = subset_mult,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = NULL
)

feature_num<-length(p.feature1)
#idents.plot <- p.feature

for (k in 1:feature_num) {
  p1 <- CoveragePlot(
    object = subset_mult,
    region = p.feature1[k],
    features = p.feature1[k],
    expression.assay = "SCT",
    extend.upstream = 500,
    extend.downstream = 10000
  )
  
  png(file=paste(genometrack_path,"/genometrack_",  p.feature[k], ".png", sep=""), width = 6, height = 5, unit="in", res=200)
  print(patchwork::wrap_plots(p1 ,ncol = 1))
  #p1
  dev.off()  
}

subset_mult
matx.atac<-na.omit(read.csv(bc_mtx))

stats_mult<-as.data.frame(cbind(nrow(subset_mult@meta.data),median(subset_mult$nCount_RNA),median(subset_mult$nFeature_RNA),median(subset_mult$nCount_ATAC),median(subset_mult$nFeature_ATAC),round(median(matx.atac$frip),digits = 2)))
colnames(stats_mult)<-c( "cell number","median UMIs per cell","median genes per cell","median fragments per cell","median number of ATAC peaks with at least one read count","median FriP")
write.csv(stats_mult,stats_multtb,quote = F)

saveRDS(subset_mult,mult_newrds)



