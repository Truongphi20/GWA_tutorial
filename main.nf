#!/usr/bin/env nextflow
include { QC_GWAS            }     from        "./src/1.QC_GWAS.nf" 
include { POP_STRATIFICATION }     from        "./src/2.Population_stratification"

workflow{
    // Inputs
    input_files_ch = channel.fromPath("./1_QC_GWAS/HapMap_3_r3_1.{bed,bim,fam}").collect()

    // Workflow
    QC_GWAS(input_files_ch)
    POP_STRATIFICATION()
}