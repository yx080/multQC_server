args <- commandArgs( trailingOnly = TRUE )

genome=args[1]
p.atac_fragments=as.numeric(args[2])
p.frip=as.numeric(args[3])
p.nfeatureRNA_max=as.numeric(args[4])
p.nfeatureRNA_min=as.numeric(args[5])
p.percent.mt=as.numeric(args[6])
p.nCountRNA=as.numeric(args[7])
i.matrix=args[8]
i.frag=args[9]
i.perbcmtx=args[10]
#sampleName=args[11]
#p.feature=args[12]
#p.feature=c(args[11],args[12],args[13],args[14],args[15],args[16])
#p.feature=c("Psca", "Krt5", "Ar", "Cd14", "Krt8", "Cd44")
o.rna=args[11]
o.atac=args[12]
o.mult=args[13]
o.bc.mtx=args[14]
o.joint_cell_calling_10x=args[15]
o.joint_cell_calling_custom=args[16]
o.venn_diagram=args[17]

#sampleName<-"f9"
#genome="hg38"
#p.atac_fragments=5000
#p.frip=0.25
#p.nfeatureRNA_max=3000000
#p.nfeatureRNA_min=20
#p.percent.mt=10
#p.nCountRNA=100
#i.matrix="/Users/yingtianxie/Desktop/testData/f9filtered_feature_bc_matrix.h5"
#i.frag="/Users/yingtianxie/Desktop/testData/f9atac_fragments.tsv.gz"
#i.perbcmtx="/Users/yingtianxie/Desktop/testData/f9per_barcode_metrics.csv"
#o.rna="analysis/f9/rna_QC/f9.scRNA.rds"
#o.atac="analysis/f9/atac_QC/f9.scATAC.rds"
#o.mult="analysis/f9/mult_QC/f9.mult.rds"
#o.bc.mtx="analysis/f9/joint_cell_calling/f9.barcodes.csv"
#o.joint_cell_calling_10x="analysis/f9/joint_cell_calling/f9.joint_cell_calling_10x.png"
#o.joint_cell_calling_custom="analysis/f9/joint_cell_calling/f9.joint_cell_calling_custom.png"
#o.venn_diagram="analysis/f9/joint_cell_calling/f9.venn_diagram.png"
#dir.create(sampleName)

libs<-""
ifelse(genome=="mm10", {libs[1] <- "EnsDb.Mmusculus.v79" ; libs[2]<-"BSgenome.Mmusculus.UCSC.mm10";},
       {libs[1]<-"EnsDb.Hsapiens.v86" ; libs[2]<-"BSgenome.Hsapiens.UCSC.hg38";})

for (i in libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}

library(ggplot2)
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(patchwork)
library(dplyr)
library(ggVennDiagram)
library(tidyr)
set.seed(1234)

counts <- Read10X_h5(i.matrix)
fragpath <- i.frag

if(genome == "mm10"){
  annotation<-GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
} else {  
  annotation<-GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
}

#head(annotation)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA data
multdata <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
multdata[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
multdata

matx<-read.csv(i.perbcmtx)
matx$umis<-matx$gex_exonic_umis + matx$gex_intronic_umis
matx$frip<-matx$atac_peak_region_fragments/matx$atac_fragments

mult.rna <- CreateSeuratObject(counts = counts$`Gene Expression`, project = "mult.rna", min.cells = 0, min.features = 0)
mult.rna[["percent.mt"]] <- PercentageFeatureSet(mult.rna, pattern = "^MT-")
mult.rna.new <- subset(mult.rna, subset = nFeature_RNA > p.nfeatureRNA_min & nFeature_RNA < p.nfeatureRNA_max & percent.mt < p.percent.mt & nCount_RNA > p.nCountRNA)
bc.rna<-mult.rna.new@meta.data

matx.atac<-subset(matx, frip > p.frip & atac_fragments > p.atac_fragments & is_cell==1)
matx.atac$cells<-"atac_cells"
matx_new<-left_join(matx,matx.atac)

matx_new$iscell<-ifelse(matx_new$is_cell==1,"yes","no")


venn_df<-list(high_quality_atac_barcodes=matx.atac$barcode,high_quality_rna_barcodes=rownames(bc.rna))
venn <- Venn(venn_df)
data <- process_data(venn)

rna_cells<-as.data.frame(rownames(bc.rna))
rna_cells$rnacells<-"rna_cells"
colnames(rna_cells)<-c("barcode","rna_cells")
matx_new2<-left_join(matx_new,rna_cells)

matx_new$atac_cells<-ifelse(matx_new$cells=='atac_cells',"yes","no")
matx_new$atac_cells<-matx_new$atac_cells %>% replace_na('no')

matx_new2$barcodes_annotation<-ifelse(matx_new2$cells=="atac_cells" & matx_new2$rna_cells=="rna_cells","cells_in_both","non_cells")
matx_new2$barcodes_annotation<-ifelse(matx_new2$cells=="atac_cells" & matx_new2$barcodes_annotation%in% NA,"atac_cells",matx_new2$barcodes_annotation)
matx_new2$barcodes_annotation<-ifelse(matx_new2$rna_cells=="rna_cells" & matx_new2$barcodes_annotation%in% NA,"rna_cells",matx_new2$barcodes_annotation)

matx_new2$iscell<-ifelse(matx_new2$is_cell==1,"yes","no")

#dir.create(paste(sampleName,"/joint_cell_calling",sep = ""))

p_cellcalling10x<-ggplot(matx_new2, aes(x=atac_peak_region_fragments, y=umis,colour=iscell)) +
  geom_point(size=0.1,alpha=0.1) + xlab("ATAC transposition events in peaks per barcode") + ylab("RNA UMIs per barcode")

#png(file=paste(sampleName,"/joint_cell_calling/", "joint_cell_calling_10x.png", sep=""), width = 3, height = 1.8, unit="in", res=200)
png(file=o.joint_cell_calling_10x, width = 3, height = 1.8, unit="in", res=200)
p_cellcalling10x+
  scale_x_log10() +
  scale_y_log10() + theme_minimal(base_size = 5)
dev.off()

p_cellcallingcustom<-ggplot(matx_new2, aes(x=atac_peak_region_fragments, y=umis,colour=barcodes_annotation)) +
  geom_point(size=0.01,alpha=0.3) + xlab("ATAC transposition events in peaks per barcode") + ylab("RNA UMIs per barcode")

#png(file=paste(sampleName,"/joint_cell_calling/", "joint_cell_calling_custom.png", sep=""), width = 3.2, height = 1.8, unit="in", res=200)
png(file=o.joint_cell_calling_custom, width = 3.2, height = 1.8, unit="in", res=200)
p_cellcallingcustom +
  scale_x_log10() +
  scale_y_log10() + theme_minimal(base_size = 5)
dev.off()

#png(file=paste(sampleName,"/joint_cell_calling/", "venn_diagram.png", sep=""), width = 3, height = 2.4, unit="in", res=200)
png(file=o.venn_diagram, width = 3, height = 2.4, unit="in", res=200)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data)) +
  scale_color_manual(values=c( "#E69F00", "#56B4E9")) +
  geom_sf(aes(color = id), data = venn_setedge(data), show.legend = FALSE) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data),size=2) +
  # 4. region label layer
  geom_sf_label(aes(label = count), data = venn_region(data),size=2) +
  theme_void() +  theme(legend.position = "none")
dev.off()

write.csv(matx_new2,file=o.bc.mtx, quote = F)

#dir.create(paste(sampleName,"/atac_QC",sep=""))
multdata1<-multdata
DefaultAssay(multdata1) <- "ATAC"
cells.use <- matx.atac$barcode
subset_atac <- subset(multdata1, cells = cells.use)
saveRDS(subset_atac, file=o.atac)

#dir.create(paste(sampleName,"/rna_QC",sep=""))
cells.use <- rownames(bc.rna)
subset_rna <- subset(mult.rna, cells = cells.use)
subset_rna
saveRDS(subset_rna, file=o.rna)

#dir.create(paste(sampleName,"/mult",sep = ""))
tb.mult<-dplyr::filter(matx_new2,barcodes_annotation=="cells_in_both")
cells.use<-tb.mult$barcode
subset_mult<- subset(multdata, cells = cells.use)
subset_mult
saveRDS(subset_mult, file=o.mult)


