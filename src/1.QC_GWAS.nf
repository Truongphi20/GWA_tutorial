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



workflow QC_GWAS {
    take:
    input_files_ch

    main:
    // Missingness
    INVESTIGATE_MISSINGNESS(input_files_ch)

    visualize_missing_script = channel.fromPath("${projectDir}/1_QC_GWAS/hist_miss.R")
    VISUALIZE_MISSINGNESS(
        INVESTIGATE_MISSINGNESS.out.output, 
        visualize_missing_script
    )
}