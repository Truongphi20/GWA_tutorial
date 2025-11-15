#!/usr/bin/env nextflow
include { QC_GWAS            }     from        "./src/1.QC_GWAS.nf" 
include { POP_STRATIFICATION }     from        "./src/2.Population_stratification"

workflow{
    // Inputs
    input_files_ch = channel.fromPath("./1_QC_GWAS/HapMap_3_r3_1.{bed,bim,fam}").collect()
    okgp_vcf = channel.fromPath("/genome_data/ALL.2of4intersection.20100804.genotypes.vcf.gz")

    // Workflow
    QC_GWAS(input_files_ch)
    POP_STRATIFICATION(
        okgp_vcf,
        QC_GWAS.out.general_qc_out,
        QC_GWAS.out.het_prune_check
    )
}