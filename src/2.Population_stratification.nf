process OKGP_CONVERT_BFILE {
    container "biocontainer/plink2:alpha2.3_jan2020"

    input:
    path okgp_vcf

    output:
    path("ALL.2of4intersection.20100804.genotypes.{bed,bim,fam}")

    script:
    """
    plink2 --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz \\
           --make-bed \\
           --max-alleles 2 \\
           --out ALL.2of4intersection.20100804.genotypes
    """
}

process OKGP_FILL_MISSING_RSID {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path original_okgp

    output:
    path(" ALL.2of4intersection.20100804.genotypes_no_missing_IDs.{bed,bim,fam}")

    script:
    """
    /usr/lib/debian-med/bin/plink --bfile ALL.2of4intersection.20100804.genotypes \\
                                  --set-missing-var-ids @:#[b37]\$1,\$2 \\
                                  --make-bed \\
                                  --out ALL.2of4intersection.20100804.genotypes_no_missing_IDs
    """
}


workflow POP_STRATIFICATION {
    take:
    okgp_vcf

    main:
    // Step 1: Prepare 1KGP data 
    OKGP_CONVERT_BFILE(okgp_vcf)
    OKGP_FILL_MISSING_RSID(OKGP_CONVERT_BFILE.out)
}