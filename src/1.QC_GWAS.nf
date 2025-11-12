process INVESTIGATE_MISSINGNESS {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path input_files

    output:
    path("plink.{lmiss,imiss}"), emit: output
    path("plink.log"), emit: logfile

    script:
    """
    export PATH=/usr/lib/debian-med/bin:/usr/lib/plink:\$PATH

    plink --bfile HapMap_3_r3_1 --missing
    """
}

process VISUALIZE_MISSINGNESS {
    container "rocker/r-base:4.5.2"

    input:
    path missingness_out
    path visualize_script

    output:
    path("histlmiss.pdf"), emit: snp_missing
    path("histimiss.pdf"), emit: sample_missing 

    script:
    """
    Rscript --no-save $visualize_script
    """
}

process FILTERING {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path input_files

    output:
    path("HapMap_3_r3_2.{bed,bim,fam}"), emit: first_filter_snp
    path("HapMap_3_r3_3.{bed,bim,fam}"), emit: first_filter_individual
    path("HapMap_3_r3_4.{bed,bim,fam}"), emit: second_filter_snp
    path("HapMap_3_r3_5.{bed,bim,fam}"), emit: second_filter_individual

    script:
    """
    # Delete SNPs with missingness >0.2
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_1 \\
                                  --geno 0.2 \\
                                  --make-bed \\
                                  --out HapMap_3_r3_2
    
    # Delete individuals with missingness >0.2.
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_2 \\
                                  --mind 0.2 \\
                                  --make-bed \\
                                  --out HapMap_3_r3_3

    # Delete SNPs with missingness >0.02.
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_3 \\
                                  --geno 0.02 \\
                                  --make-bed \\
                                  --out HapMap_3_r3_4
    
    # Delete individuals with missingness >0.02.
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_4 \\
                                  --mind 0.02 \\
                                  --make-bed \\
                                  --out HapMap_3_r3_5
    """
}



workflow QC_GWAS {
    take:
    input_files_ch     // file: [./1_QC_GWAS/HapMap_3_r3_1.{bed,bim,fam}]

    main:
    // Survey Missingness
    INVESTIGATE_MISSINGNESS(input_files_ch)

    visualize_missing_script = channel.fromPath("${projectDir}/1_QC_GWAS/hist_miss.R")
    VISUALIZE_MISSINGNESS(
        INVESTIGATE_MISSINGNESS.out.output, 
        visualize_missing_script
    )

    // Filtering
    FILTERING(input_files_ch)


}