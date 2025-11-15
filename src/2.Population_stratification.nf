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
    container "biocontainer/plink2:alpha2.3_jan2020"

    input:
    path original_okgp

    output:
    path("ALL.2of4intersection.20100804.genotypes_no_missing_IDs.{bed,bim,fam}")

    script:
    """
    plink2 --bfile ALL.2of4intersection.20100804.genotypes \\
            --set-missing-var-ids @:#[b37]\\\$r,\\\$a \\
            --make-bed \\
            --out ALL.2of4intersection.20100804.genotypes_no_missing_IDs
    """
}

process OKGP_QC_DROP_VARIANTS {
    container "biocontainer/plink2:alpha2.3_jan2020"

    input:
    path no_missing_ids

    output:
    path("1kG_MDS.{bed,bim,fam}") 

    script:
    """
    plink2 --bfile ALL.2of4intersection.20100804.genotypes_no_missing_IDs \\
          --geno 0.2 \\
          --allow-no-sex \\
          --make-bed \\
          --out 1kG_MDS
    """
}

process OKGP_QC_DROP_INDIVIDUALS {
    container "biocontainer/plink2:alpha2.3_jan2020"

    input:
    path drop_variants

    output:
    path("1kG_MDS2.{bed,bim,fam}") 

    script:
    """
    plink2 --bfile 1kG_MDS \\
           --mind 0.2 \\
           --allow-no-sex \\
           --make-bed \\
           --out 1kG_MDS2
    """
}

process OKGP_QC_STRICTHEN_DROP_VARIANTS {
    container "biocontainer/plink2:alpha2.3_jan2020"

    input:
    path no_missing_ids

    output:
    path("1kG_MDS3.{bed,bim,fam}") 

    script:
    """
    plink2 --bfile 1kG_MDS2 \\
          --geno 0.02 \\
          --allow-no-sex \\
          --make-bed \\
          --out 1kG_MDS3
    """
}

process OKGP_QC_STRICTHEN_DROP_INDIVIDUALS {
    container "biocontainer/plink2:alpha2.3_jan2020"

    input:
    path drop_variants

    output:
    path("1kG_MDS4.{bed,bim,fam}") 

    script:
    """
    plink2 --bfile 1kG_MDS3 \\
           --mind 0.02 \\
           --allow-no-sex \\
           --make-bed \\
           --out 1kG_MDS4
    """
}

process OKGP_QC_MAF {
    container "biocontainer/plink2:alpha2.3_jan2020"

    input:
    path drop_vars_and_ind

    output:
    path("1kG_MDS5.{bed,bim,fam}")

    script:
    """
    plink2 --bfile 1kG_MDS4 \\
            --maf 0.05 \\
            --allow-no-sex \\
            --make-bed \\
            --out 1kG_MDS5
    """
}


workflow POP_STRATIFICATION {
    take:
    okgp_vcf

    main:
    // Step 1: Prepare 1KGP data 
    OKGP_CONVERT_BFILE(okgp_vcf)
    OKGP_FILL_MISSING_RSID(OKGP_CONVERT_BFILE.out)
    OKGP_QC_DROP_VARIANTS(OKGP_FILL_MISSING_RSID.out)
    OKGP_QC_DROP_INDIVIDUALS(OKGP_QC_DROP_VARIANTS.out)
    OKGP_QC_STRICTHEN_DROP_VARIANTS(OKGP_QC_DROP_INDIVIDUALS.out)
    OKGP_QC_STRICTHEN_DROP_INDIVIDUALS(OKGP_QC_STRICTHEN_DROP_VARIANTS.out)
    OKGP_QC_MAF(OKGP_QC_STRICTHEN_DROP_INDIVIDUALS.out)
}