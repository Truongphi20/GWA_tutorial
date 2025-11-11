#!/usr/bin/env nextflow
include { QC_GWAS }     from        "./src/1.QC_GWAS.nf" 

workflow{
    QC_GWAS()
}