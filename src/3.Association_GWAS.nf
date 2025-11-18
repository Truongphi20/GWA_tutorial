process ANALYSIS_ASSOC {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path associate_bfile

    output:
    path("assoc_results.assoc")

    script:
    """
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_13 --assoc --out assoc_results
    """
}

process ANALYSIS_LOGISTIC {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path associate_bfile
    path covar_mds

    output:
    path("logistic_results.assoc.logistic")

    script:
    """
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_13 \\
                                  --covar covar_mds.txt \\
                                  --logistic \\
                                  --hide-covar \\
                                  --out logistic_results
    """
}


workflow ASSOCIATION_GWAS {
    take:
    associate_bfile                 // HapMap_3_r3_13.{bed,bim,fam}
    covar_mds                       // covar_mds.txt

    main:
    ANALYSIS_ASSOC(associate_bfile)
    ANALYSIS_LOGISTIC(
        associate_bfile,
        covar_mds
    )
    
}