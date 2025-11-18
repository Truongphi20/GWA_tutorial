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

process TESTING_ADJUSTING {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path associate_bfile

    output:
    path("adjusted_assoc_results.assoc")

    script:
    """
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_13 \\
                                  --assoc --adjust \\
                                  --out adjusted_assoc_results
    """
}

process TESTING_PERMUTATION {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path associate_bfile

    output:
    path("sorted_subset.txt")

    script:
    """
    # Generate subset of SNPs
    awk '{ if (\$4 >= 21595000 && \$4 <= 21605000) print \$2 }' HapMap_3_r3_13.bim > subset_snp_chr_22.txt
    
    # Filter your bfile based on the subset of SNPs generated in the step above.
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_13 \\
                                  --extract subset_snp_chr_22.txt \\
                                  --make-bed \\
                                  --out HapMap_subset_for_perm

    # Perform 1000000 perrmutations.
    /usr/lib/debian-med/bin/plink --bfile HapMap_subset_for_perm \\
                                  --assoc \\
                                  --mperm 1000000 \\
                                  --out subset_1M_perm_result

    # Order your data, from lowest to highest p-value.
    sort -gk 4 subset_1M_perm_result.assoc.mperm > sorted_subset.txt
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

    // Step 2: Testing
    TESTING_ADJUSTING(associate_bfile)
    TESTING_PERMUTATION(associate_bfile)

    
}