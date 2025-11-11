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



workflow QC_GWAS {
    take:
    input_files_ch

    main:
    INVESTIGATE_MISSINGNESS(input_files_ch)
}