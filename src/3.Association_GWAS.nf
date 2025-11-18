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

process ANALYSIS_CLEANUP {

    input:
    path logistic_output

    output:
    path("logistic_results.assoc_2.logistic")

    script:
    """
    awk '!/'NA'/' logistic_results.assoc.logistic > logistic_results.assoc_2.logistic
    """
}


workflow ASSOCIATION_GWAS {
    take:
    associate_bfile                 // HapMap_3_r3_13.{bed,bim,fam}
    covar_mds                       // covar_mds.txt

    main:
    // Step 1: Associate analysis 
    ANALYSIS_ASSOC(associate_bfile)
    ANALYSIS_LOGISTIC(
        associate_bfile,
        covar_mds
    )
    ANALYSIS_CLEANUP(ANALYSIS_LOGISTIC.out)

    
}