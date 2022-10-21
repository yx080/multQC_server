configfile: "config.yaml"

# Check metasheet setup
from cmath import log
from random import sample
from scripts.metasheet_setup import updateMeta
import pandas as pd
import yaml
import os


config = updateMeta(config)
metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#')

pwd = os.getcwd()

rule download_example_GR_ChIP:
    params:
        hm_data = "https://www.dropbox.com/s/53hj7hv02pcxde0/testhumandata.tar.gz?dl=0",
    output:
        "test_data/hm/"
    shell:
        "wget {params.hm_data} -O hmtest_data.tar.gz && tar -xvzf hmtest_data.tar.gz "

rule make_matrix_folder:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "data/matrix/{sample}.h5"
    shell:
        "cp {input} {output}"

rule make_fragmentfile_folder:
    input:
        lambda wildcards: config["frag"][wildcards.sample]
    output:
        "data/frag/{sample}.tsv.gz"
    shell:
        "cp {input} {output}"

rule make_fragmentfileindex_folder:
    input:
        lambda wildcards: config["frag_index"][wildcards.sample]
    output:
        "data/frag/{sample}.tsv.gz.tbi"
    shell:
        "cp {input} {output}"

rule make_perbcmtx_folder:
    input:
        lambda wildcards: config["perbcmtx"][wildcards.sample]
    output:
        "data/perbcmtx/{sample}.csv"
    shell:
        "cp {input} {output}"  

rule make_summary_folder:
    input:
        lambda wildcards: config["summary"][wildcards.sample]
    output:
        "data/summary/{sample}.csv"
    shell:
        "cp {input} {output}"  

rule generate_obj:
    input:
        #matrix=expand("data/matrix/{sample}.h5", sample=config["samples"]),
        #frag=expand("data/frag/{sample}.tsv.gz", sample=config["samples"]),
        #perbcmtx=expand("data/perbcmtx/{sample}.csv", sample=config["samples"]),
        #frag_index=expand("data/frag/{sample}.tsv.gz.tbi", sample=config["samples"]),
        matrix="data/matrix/{sample}.h5",
        frag="data/frag/{sample}.tsv.gz",
        perbcmtx="data/perbcmtx/{sample}.csv",
        frag_index="data/frag/{sample}.tsv.gz.tbi",
    params:
        #dataType=config["dataType"],
        genome=config["assembly"],
        atac_fragments=config["atac_fragments"],
        frip=config["frip"],
        nfeatureRNA_max=config["nfeatureRNA_max"],
        nfeatureRNA_min=config["nfeatureRNA_min"],
        percent_mt=config["percent_mt"],
        nCountRNA=config["nCountRNA"],
    output: 
        rnards="analysis/{sample}/rna_QC/{sample}.scRNA.rds",
        atacrds="analysis/{sample}/atac_QC/{sample}.scATAC.rds",
        multrds="analysis/{sample}/mult_QC/{sample}.mult.rds",
        bc_mtx="analysis/{sample}/joint_cell_calling/{sample}.barcodes.csv",
        joint_cell_calling_10x="analysis/{sample}/joint_cell_calling/{sample}.joint_cell_calling_10x.png",
        joint_cell_calling_custom="analysis/{sample}/joint_cell_calling/{sample}.joint_cell_calling_custom.png",
        venn_diagram="analysis/{sample}/joint_cell_calling/{sample}.venn_diagram.png"
    log:
        "log/generate_obj.{sample}.log"
    shell:
        "Rscript scripts/filter_obj.R {params.genome} {params.atac_fragments} {params.frip} {params.nfeatureRNA_max} {params.nfeatureRNA_min} {params.percent_mt} {params.nCountRNA} "
        "{input.matrix} {input.frag} {input.perbcmtx} {output.rnards} {output.atacrds} {output.multrds} {output.bc_mtx} {output.joint_cell_calling_10x} "
        "{output.joint_cell_calling_custom} {output.venn_diagram}"
        #"""
        #if [ "{params.dataType}" == "multiomic" ]; then  (Rscript generate_results.new.R {params.genome} {params.atac_fragments} {params.frip} {params.nfeatureRNA_max} {params.nfeatureRNA_min} {params.percent_mt} {params.nCountRNA} {input.matrix} {input.frag} {input.perbcmtx} {params.sample} ) >{output}
        #elif [ "{params.dataType}" == "scATAC" ]; then (Rscript generate_scatac_results.R {params.genome} {params.atac_fragments} {params.frip} {params.nfeatureRNA_max} {params.nfeatureRNA_min} {params.percent_mt} {params.nCountRNA} {input.matrix} {input.frag} {input.perbcmtx} 'testAP' ) >{output}
        #else (Rscript generate_scrna_results.R {params.genome} {params.atac_fragments} {params.frip} {params.nfeatureRNA_max} {params.nfeatureRNA_min} {params.percent_mt} {params.nCountRNA} {input.matrix} {input.frag} {input.perbcmtx} 'testAP' ) >{output}
        #fi
        #"""


# rule scatac_qc:
#     input:
#         atacrds="analysis/{sample}/atac_QC/{sample}.scATAC.rds",
#         bc_mtx="analysis/{sample}/joint_cell_calling/{sample}.barcodes.csv",
#     params:
#         feature=config["feature"],
#         genome=config["assembly"],
#         motif_path="analysis/{sample}/atac_QC/motif_analysis",
#         diff_path="analysis/{sample}/atac_QC/diff_analysis",
#     output:
#         peak_target="analysis/{sample}/atac_QC/{sample}.peak_target.png",
#         tss_insertsize="analysis/{sample}/atac_QC/{sample}.tss_insertsize.png",
#         atac_qc_violin="analysis/{sample}/atac_QC/{sample}.atac_qc_violin.png",
#         atac_depthCor="analysis/{sample}/atac_QC/{sample}.atac_depthCor.png",
#         atac_umap_cluster="analysis/{sample}/atac_QC/{sample}.atac_umap_cluster.png",
#         atac_coordinates="analysis/{sample}/atac_QC/{sample}.atac_coordinates.csv",
#         atac_umap_byvennbarcode="analysis/{sample}/atac_QC/{sample}.atac_umap_byvennbarcode.png",
#         atac_umap_depth="analysis/{sample}/atac_QC/{sample}.atac_umap_depth.png",
#         atac_umap_marker_gene="analysis/{sample}/atac_QC/{sample}.atac_umap_marker_gene.png",          
#     log:
#         "log/scatac_qc.{sample}.log"
#     shell:
#         "Rscript scripts/scATAC_QC.R {input.atacrds} {input.bc_mtx} {output.peak_target} {output.tss_insertsize} {output.atac_qc_violin} {output.atac_depthCor} {output.atac_umap_cluster} "
#         "{output.atac_coordinates} {output.atac_umap_byvennbarcode} {output.atac_umap_depth} {params.diff_path} {output.atac_umap_marker_gene} {params.feature} {params.motif_path} "
#         "{params.genome}"  

rule scatac_qcplots_qcumap:
    input:
        atacrds="analysis/{sample}/atac_QC/{sample}.scATAC.rds",
        bc_mtx="analysis/{sample}/joint_cell_calling/{sample}.barcodes.csv",
    params:
     #   feature=config["feature"],
        genome=config["assembly"],
    #    motif_path="analysis/{sample}/atac_QC/motif_analysis",
     #   diff_path="analysis/{sample}/atac_QC/diff_analysis",
    output:
        peak_target="analysis/{sample}/atac_QC/{sample}.peak_target.png",
        tss_insertsize="analysis/{sample}/atac_QC/{sample}.tss_insertsize.png",
        atac_qc_violin="analysis/{sample}/atac_QC/{sample}.atac_qc_violin.png",
        atac_depthCor="analysis/{sample}/atac_QC/{sample}.atac_depthCor.png",
        atac_umap_cluster="analysis/{sample}/atac_QC/{sample}.atac_umap_cluster.png",
        atac_coordinates="analysis/{sample}/atac_QC/{sample}.atac_coordinates.csv",
        atac_umap_byvennbarcode="analysis/{sample}/atac_QC/{sample}.atac_umap_byvennbarcode.png",
        atac_umap_depth="analysis/{sample}/atac_QC/{sample}.atac_umap_depth.png", 
        atacqcrds="analysis/{sample}/atac_QC/{sample}.scATAC.qc.rds"       
    log:
        "log/scatac_qc.{sample}.log"
    shell:
        "Rscript scripts/scATAC_QC.R {input.atacrds} {input.bc_mtx} {output.peak_target} {output.tss_insertsize} {output.atac_qc_violin} {output.atac_depthCor} {output.atac_umap_cluster} "
        "{output.atac_coordinates} {output.atac_umap_byvennbarcode} {output.atac_umap_depth} {params.genome} {output.atacqcrds}"
 

rule scatac_diff_analysis_markergene:
    input:
        atacrds="analysis/{sample}/atac_QC/{sample}.scATAC.qc.rds",
    #    bc_mtx="analysis/{sample}/joint_cell_calling/{sample}.barcodes.csv",
    params:
        feature=config["feature"],
        genome=config["assembly"],
    #    motif_path="analysis/{sample}/atac_QC/motif_analysis",
        diff_path="analysis/{sample}/atac_QC/diff_analysis",
    output:
        atac_umap_marker_gene="analysis/{sample}/atac_QC/{sample}.atac_umap_marker_gene.png",          
    log:
        "log/scatac_df.{sample}.log"
    shell:
        "Rscript scripts/scATAC_diff_analysis_markergene.R {input.atacrds} {params.diff_path} {output.atac_umap_marker_gene} {params.feature} {params.genome}"
        

rule scatac_motif:
    input:
        atacrds="analysis/{sample}/atac_QC/{sample}.scATAC.qc.rds",
     #   bc_mtx="analysis/{sample}/joint_cell_calling/{sample}.barcodes.csv",
    params:
     #   feature=config["feature"],
        genome=config["assembly"],
        motif_path="analysis/{sample}/atac_QC/motif_analysis",
    #    diff_path="analysis/{sample}/atac_QC/diff_analysis",
    output:
        outfile1="analysis/{sample}/atac_QC/motif_analysis/motif_analysis_cluster_1.csv"              
    log:
        "log/scatac_mtf.{sample}.log"
    shell:
        "Rscript scripts/scATAC_motif.R {input.atacrds} {params.motif_path} {params.genome} {output.outfile1}"

rule all:
    input:
        expand("analysis/{sample}/mult_QC/{sample}.mult.rds",sample=config["samples"]),
        expand("analysis/{sample}/atac_QC/{sample}.scATAC.qc.rds",sample=config["samples"]),
        expand("analysis/{sample}/atac_QC/{sample}.atac_umap_marker_gene.png",sample=config["samples"]),
        expand("analysis/{sample}/atac_QC/motif_analysis/motif_analysis_cluster_1.csv",sample=config["samples"]),
        expand("analysis/{sample}/mult_QC/{sample}.stats_mult.csv",sample=config["samples"]),
        expand("analysis/{sample}/joint_cell_calling/{sample}.venn_diagram.png",sample=config["samples"]),
        expand("analysis/{sample}/rna_QC/{sample}.rna_heatmap.png",sample=config["samples"]),
        expand("analysis/{sample}/rna_QC/{sample}.rna_cell_cycle.png",sample=config["samples"]),
        expand("data/summary/{sample}.csv",sample=config["samples"]),
        expand("log/report.{sample}.log",sample=config["samples"])

# rule scrna_qc:
#     input:   
#         rnards="analysis/{sample}/rna_QC/{sample}.scRNA.rds",
#         bc_mtx="analysis/{sample}/joint_cell_calling/{sample}.barcodes.csv",
#         matrix="data/matrix/{sample}.h5",
#     params:
#         diff_path="analysis/{sample}/rna_QC/diff_analysis",
#         p_feature=config["feature"],
#     output:    
#         rna_qc_violin="analysis/{sample}/rna_QC/{sample}.rna_qc_violin.png",
#         feature_cor="analysis/{sample}/rna_QC/{sample}.feature_cor.png",
#         rna_highly_variable_gene="analysis/{sample}/rna_QC/{sample}.rna_highly_variable_gene.png",
#         rna_umap_cluster="analysis/{sample}/rna_QC/{sample}.rna_umap_cluster.png",
#         rna_coordinates="analysis/{sample}/rna_QC/{sample}.rna_coordinates.csv",
#         rna_umap_byvennbarcode="analysis/{sample}/rna_QC/{sample}.rna_umap_byvennbarcode.png",
#         rna_umap_depth="analysis/{sample}/rna_QC/{sample}.rna_umap_depth.png",
#         rna_umap_mt_percent="analysis/{sample}/rna_QC/{sample}.rna_umap_mt_percent.png",
#         rna_umap_cell_cycle="analysis/{sample}/rna_QC/{sample}.rna_umap_cell_cycle.png",
#         rna_cell_cycle="analysis/{sample}/rna_QC/{sample}.rna_cell_cycle.png",
#         diff_results="analysis/{sample}/rna_QC/diff_analysis/{sample}.diff_results.csv",
#         rna_ump_markergene="analysis/{sample}/rna_QC/{sample}.rna_ump_markergene.png",
#         rna_heatmap="analysis/{sample}/rna_QC/{sample}.rna_heatmap.png",
#     log:
#         "log/scrna.qc.{sample}.log"
#     shell:
#         "Rscript scripts/scRNA_QC.R {input.rnards} {output.rna_qc_violin} {output.feature_cor} {output.rna_highly_variable_gene} {output.rna_umap_cluster} {output.rna_coordinates} "
#         "{output.rna_umap_byvennbarcode} {output.rna_umap_depth} {output.rna_umap_mt_percent} {output.rna_umap_cell_cycle} {output.rna_cell_cycle} {output.diff_results} {params.diff_path} "
#         "{params.p_feature} {output.rna_ump_markergene} {output.rna_heatmap} {input.bc_mtx} {input.matrix}"

rule scrna_qcplots_qcumap:
    input:   
        rnards="analysis/{sample}/rna_QC/{sample}.scRNA.rds",
        bc_mtx="analysis/{sample}/joint_cell_calling/{sample}.barcodes.csv",
        matrix="data/matrix/{sample}.h5",
    #params:
        #diff_path="analysis/{sample}/rna_QC/diff_analysis",
        #p_feature=config["feature"],
    output:    
        rna_qc_violin="analysis/{sample}/rna_QC/{sample}.rna_qc_violin.png",
        feature_cor="analysis/{sample}/rna_QC/{sample}.feature_cor.png",
        rna_highly_variable_gene="analysis/{sample}/rna_QC/{sample}.rna_highly_variable_gene.png",
        rna_umap_cluster="analysis/{sample}/rna_QC/{sample}.rna_umap_cluster.png",
        rna_coordinates="analysis/{sample}/rna_QC/{sample}.rna_coordinates.csv",
        rna_umap_byvennbarcode="analysis/{sample}/rna_QC/{sample}.rna_umap_byvennbarcode.png",
        rna_umap_depth="analysis/{sample}/rna_QC/{sample}.rna_umap_depth.png",
        rna_umap_mt_percent="analysis/{sample}/rna_QC/{sample}.rna_umap_mt_percent.png",
        rnaqcrds="analysis/{sample}/rna_QC/{sample}.scRNA.qc.rds",
    log:
        "analysis/{sample}/rna_QC/{sample}.scRNA.qc.rds"
    shell:
        "Rscript scripts/scRNA_QC.R {input.rnards} {output.rna_qc_violin} {output.feature_cor} {output.rna_highly_variable_gene} {output.rna_umap_cluster} {output.rna_coordinates} "
        "{output.rna_umap_byvennbarcode} {output.rna_umap_depth} {output.rna_umap_mt_percent} {input.bc_mtx} {input.matrix} {output.rnaqcrds}"

rule scrna_diff_analysis:
    input:   
        rnards="analysis/{sample}/rna_QC/{sample}.scRNA.qc.rds",
    #    bc_mtx="analysis/{sample}/joint_cell_calling/{sample}.barcodes.csv",
    #    matrix="data/matrix/{sample}.h5",
    params:
        diff_path="analysis/{sample}/rna_QC/diff_analysis",
        p_feature=config["feature"],
    output:    
        diff_results="analysis/{sample}/rna_QC/diff_analysis/{sample}.diff_results.csv",
        rna_ump_markergene="analysis/{sample}/rna_QC/{sample}.rna_ump_markergene.png",
        rna_heatmap="analysis/{sample}/rna_QC/{sample}.rna_heatmap.png",
    log:
        "log/scrna.df.{sample}.log"
    shell:
        "Rscript scripts/scRNA_diff_analysis.R {input.rnards} {output.diff_results} {params.diff_path} {params.p_feature} {output.rna_ump_markergene} {output.rna_heatmap} "

rule scrna_qc_cell_cycle:
    input:   
        rnards="analysis/{sample}/rna_QC/{sample}.scRNA.qc.rds",
        bc_mtx="analysis/{sample}/joint_cell_calling/{sample}.barcodes.csv",
        matrix="data/matrix/{sample}.h5",
    #params:
    #    diff_path="analysis/{sample}/rna_QC/diff_analysis",
    #    p_feature=config["feature"],
    output:    
        rna_umap_cell_cycle="analysis/{sample}/rna_QC/{sample}.rna_umap_cell_cycle.png",
        rna_cell_cycle="analysis/{sample}/rna_QC/{sample}.rna_cell_cycle.png",
    log:
        "log/scrna.cc.{sample}.log"
    shell:
        "Rscript scripts/scRNA_cellcycle.R {input.rnards} {output.rna_umap_cell_cycle} {output.rna_cell_cycle}  "
        "{input.bc_mtx} {input.matrix}"

rule mult_qc:
    input:
        multrds="analysis/{sample}/mult_QC/{sample}.mult.rds",
        bc_mtx="analysis/{sample}/joint_cell_calling/{sample}.barcodes.csv",
    params:
        genome=config["assembly"],
        p_feature=config["feature"],
        genometrack_path="analysis/{sample}/mult_QC/",
    output:
        jointumap="analysis/{sample}/mult_QC/{sample}.jointumap.png",
        mult_coordinates="analysis/{sample}/mult_QC/{sample}.mult_coordinates.csv",
        jointumap_rna="analysis/{sample}/mult_QC/{sample}.jointumap_rna.png",
        jointumap_atac="analysis/{sample}/mult_QC/{sample}.jointumap_atac.png",
     #   stats_mult="analysis/{sample}/mult_QC/{sample}.stats_mult.csv",
        mult_outrds="analysis/{sample}/mult_QC/{sample}.out.mult.rds"
    log:
        "log/mult.qc.{sample}.log"
    shell:
        "Rscript scripts/mult_QC.orig.R {input.multrds} {output.jointumap} {output.mult_coordinates} {output.jointumap_rna} {output.jointumap_atac} {params.genometrack_path} "
        "{params.genome} {params.p_feature} {input.bc_mtx} {output.mult_outrds}"

rule mult_coverageplot:
    input:
        multoutrds="analysis/{sample}/mult_QC/{sample}.out.mult.rds",
        bc_mtxcov="analysis/{sample}/joint_cell_calling/{sample}.barcodes.csv",
    params:
        genome=config["assembly"],
        p_feature=config["feature"],
        genometrack_path="analysis/{sample}/mult_QC/",
    output:
        stats_mult="analysis/{sample}/mult_QC/{sample}.stats_mult.csv",
        mult_outrds="analysis/{sample}/mult_QC/{sample}.cov.out.mult.rds",
    log:
        "log/mult.cov.{sample}.log"
    shell:
        "Rscript scripts/mult_coverageplot.R {input.multoutrds} {params.genometrack_path} {output.stats_mult} {params.genome} {params.p_feature} {output.mult_outrds} "
        "{input.bc_mtxcov}"    


rule generate_reports:
    input:
        summary_tb="data/summary/{sample}.csv",
        stats_mult_tb="analysis/{sample}/mult_QC/{sample}.stats_mult.csv",    
    params:     
        path= pwd,
        sampleName="{sample}",
        atac_fragments=config["atac_fragments"],
        frip=config["frip"],
        nfeatureRNA_max=config["nfeatureRNA_max"],
        nfeatureRNA_min=config["nfeatureRNA_min"],
        percent_mt=config["percent_mt"],
        nCountRNA=config["nCountRNA"],
        p_feature=config["feature"],
    #output:
        #report="analysis/reports/{sample}.reports.html"
    log:
        "log/report.{sample}.log"
    shell:
        #"export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/quarto/bin && "
        "Rscript -e \"rmarkdown::render('scripts/generate_report.bak.Rmd',params=list(summary_tb=\'{input.summary_tb}\', stats_mult_tb=\'{input.stats_mult_tb}\', sampleName=\'{params.sampleName}\', pwd=\'{params.path}\', patac_fragments={params.atac_fragments}, pfrip={params.frip}, pnfeatureRNA_max={params.nfeatureRNA_max}, pnfeatureRNA_min={params.nfeatureRNA_min}, ppercentmt={params.percent_mt}, pnCountRNA={params.nCountRNA},pfeature=\'{params.p_feature}\'),output_file='{params.path}/analysis/{params.sampleName}/reports.html')\" "
        #"Rscript -e \"rmarkdown::render('scripts/generate_report.bak.Rmd')\" "
        
