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

process HAMONIZE_OKGP_VARIANT {
    container "biocontainer/plink2:alpha2.3_jan2020"

    input:
    path general_qc_out
    path drop_by_maf

    output:
    path("1kG_MDS6.{bed,bim,fam}")

    script:
    """
    awk '{print\$2}' HapMap_3_r3_12.bim > HapMap_SNPs.txt
    plink2 --bfile 1kG_MDS5 \\
            --extract HapMap_SNPs.txt \\
            --make-bed \\
            --out 1kG_MDS6
    """
}

process HAMONIZE_HAPMAP_VARIANT {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path hamonized_okgp
    path general_qc_out

    output:
    path("HapMap_MDS.{bed,bim,fam}")

    script:
    """
    awk '{print\$2}' 1kG_MDS6.bim > 1kG_MDS6_SNPs.txt
    /usr/lib/debian-med/bin/plink \\
           --bfile HapMap_3_r3_12 \\
           --extract 1kG_MDS6_SNPs.txt \\
           --recode \\
           --make-bed \\
           --out HapMap_MDS
    """
}

process HAMONIZE_OKGP_BUILD {
    container "biocontainer/plink2:alpha2.3_jan2020"

    input:
    path okgp_harmonize_var
    path hapmap_harmonize_var

    output:
    path("1kG_MDS7.{bed,bim,fam}")

    script:
    """
    awk '{print\$2,\$4}' HapMap_MDS.bim > buildhapmap.txt
    plink2 --bfile 1kG_MDS6 \\
           --update-map buildhapmap.txt \\
           --make-bed \\
           --out 1kG_MDS7
    """
}

process MERGE_ENSURE_HAPMAP_REF {
    container "biocontainer/plink2:alpha2.3_jan2020"

    input:
    path okgp_harmonize_build
    path hapmap_harmonize_var

    output:
    path("HapMap-adj.{bed,bim,fam}"), emit: bfiles
    path("1kg_ref-list.txt"), emit: ref_list

    script:
    """
    awk '{print\$2,\$5}' 1kG_MDS7.bim > 1kg_ref-list.txt
    plink2 --bfile HapMap_MDS \\
           --ref-allele 1kg_ref-list.txt \\
           --make-bed \\
           --out HapMap-adj
    """
}

process MERGE_STRAND_PB_CHECK {

    input:
    path okgp_harmonize_build
    path hapmap_ensure_ref

    output:
    path("all_differences.txt")

    script:
    """
    awk '{print\$2,\$5,\$6}' 1kG_MDS7.bim > 1kGMDS7_tmp
    awk '{print\$2,\$5,\$6}' HapMap-adj.bim > HapMap-adj_tmp
    sort 1kGMDS7_tmp HapMap-adj_tmp |uniq -u > all_differences.txt
    """
} 

process MERGE_STRAND_PB_FLIP {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path strand_check
    path hapmap_ensure_ref
    path hapmap_ensure_ref_list

    output:
    path("corrected_hapmap.{bed,bim,fam}"), emit: bfiles 
    path("uncorresponding_SNPs.txt"), emit: promblematic_still

    script:
    """
    awk '{print\$1}' all_differences.txt | sort -u > flip_list.txt
    /usr/lib/debian-med/bin/plink \\
            --bfile HapMap-adj \\
            --flip flip_list.txt \\
            --reference-allele 1kg_ref-list.txt \\
            --make-bed \\
            --out corrected_hapmap

    awk '{print\$2,\$5,\$6}' corrected_hapmap.bim > corrected_hapmap_tmp
    sort 1kGMDS7_tmp corrected_hapmap_tmp |uniq -u  > uncorresponding_SNPs.txt
    """
}

process MERGE_STRAND_PB_DROP_VARIANTS {
    container "biocontainer/plink2:alpha2.3_jan2020"

    input:
    path corrected_hapmap
    path okgp_harmonize_build

    output:
    path("HapMap_MDS2.{bed,bim,fam}")       , emit: hapmap 
    path("1kG_MDS8.{bed,bim,fam}")          , emit: okgp

    script:
    """
    awk '{print\$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
    plink2 --bfile corrected_hapmap \\
           --exclude SNPs_for_exlusion.txt \\
           --make-bed \\
           --out HapMap_MDS2
    
    plink2 --bfile 1kG_MDS7 \\
           --exclude SNPs_for_exlusion.txt \\
           --make-bed \\
           --out 1kG_MDS8
    """
}


workflow POP_STRATIFICATION {
    take:
    okgp_vcf
    general_qc_out              // HapMap_3_r3_12.fam,HapMap_3_r3_12.bed, and HapMap_3_r3_12.bim
    het_prune_check             // indepSNP.prune.in

    main:
    // Step 1: Prepare 1KGP data 
    OKGP_CONVERT_BFILE(okgp_vcf)
    OKGP_FILL_MISSING_RSID(OKGP_CONVERT_BFILE.out)
    OKGP_QC_DROP_VARIANTS(OKGP_FILL_MISSING_RSID.out)
    OKGP_QC_DROP_INDIVIDUALS(OKGP_QC_DROP_VARIANTS.out)
    OKGP_QC_STRICTHEN_DROP_VARIANTS(OKGP_QC_DROP_INDIVIDUALS.out)
    OKGP_QC_STRICTHEN_DROP_INDIVIDUALS(OKGP_QC_STRICTHEN_DROP_VARIANTS.out)
    OKGP_QC_MAF(OKGP_QC_STRICTHEN_DROP_INDIVIDUALS.out)

    // Step 2: Harmonize 1KGP and HapMap data
    HAMONIZE_OKGP_VARIANT(
        general_qc_out,
        OKGP_QC_MAF.out
    )

    HAMONIZE_HAPMAP_VARIANT(
        HAMONIZE_OKGP_VARIANT.out,
        general_qc_out
    )

    HAMONIZE_OKGP_BUILD(
        HAMONIZE_OKGP_VARIANT.out,
        HAMONIZE_HAPMAP_VARIANT.out
    )

    // Step 3: Merge 1KGP and HapMap data
    MERGE_ENSURE_HAPMAP_REF(
        HAMONIZE_OKGP_BUILD.out,
        HAMONIZE_HAPMAP_VARIANT.out
    )

    MERGE_STRAND_PB_CHECK(
        HAMONIZE_OKGP_BUILD.out,
        MERGE_ENSURE_HAPMAP_REF.out.bfiles
    )

    MERGE_STRAND_PB_FLIP(
        MERGE_STRAND_PB_CHECK.out,
        MERGE_ENSURE_HAPMAP_REF.out.bfiles,
        MERGE_ENSURE_HAPMAP_REF.out.ref_list
    )

    MERGE_STRAND_PB_DROP_VARIANTS(
        MERGE_STRAND_PB_FLIP.out.bfiles,
        HAMONIZE_OKGP_BUILD.out
    )

}