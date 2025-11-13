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

process FILTERING_MISSINGNESS {
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

process CHECK_SEX_DISCREPANCY {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path filter_missing_output

    output:
    path("plink.sexcheck")

    script:
    """
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_5 --check-sex
    """
}

process SEX_DISCREPANCY_VISUALIZE {
    container "rocker/r-base:4.5.2"

    input:
    path sex_check
    path visualize_script

    output:
    path("Gender_check.pdf"), emit: gender
    path("Men_check.pdf"), emit: men
    path("Women_check.pdf"), emit: women

    script:
    """
    Rscript --no-save $visualize_script
    """
}

process SEX_DISCREPANCY_HANDLE_BY_REMOVE {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path sex_check
    path filter_missing_output

    output:
    path("HapMap_3_r3_6.{bed,bim,fam}")

    script:
    """
    grep "PROBLEM" $sex_check| awk '{print\$1,\$2}'> sex_discrepancy.txt
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_5  \\
                                  --remove sex_discrepancy.txt \\
                                  --make-bed \\
                                  --out HapMap_3_r3_6 
    """
}

process SEX_DISCREPANCY_HANDLE_BY_IMPUTE {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path filter_missing_output

    output:
    path("HapMap_3_r3_6.{bed,bim,fam}")

    script:
    """
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_5 \\
                                    --impute-sex \\
                                    --make-bed \\
                                    --out HapMap_3_r3_6
    """
}

process MAF_AUTOSOMAL_SELECTION {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path sex_discrepancy_output

    output:
    path("HapMap_3_r3_7.{bed,bim,fam}"), emit: plink_files
    path("snp_1_22.txt")               , emit: snp_list
    path("MAF_check.frq")              , emit: check_frq

    script:
    """
    awk '{ if (\$1 >= 1 && \$1 <= 22) print \$2 }' HapMap_3_r3_6.bim > snp_1_22.txt
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_6 \\
                                  --extract snp_1_22.txt \\
                                  --make-bed \\
                                  --out HapMap_3_r3_7

    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_7 \\
                                  --freq \\
                                  --out MAF_check
    """
}

process MAF_PLOT_DISTRIBUTION {
    container "rocker/r-base:4.5.2"

    input:
    path check_frq_maf
    path maf_check_script

    output:
    path("MAF_distribution.pdf")

    script:
    """
    Rscript --no-save $maf_check_script
    """

}

process MAF_FILTERING {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path maf_selection_output

    output:
    path("HapMap_3_r3_8.{bed,bim,fam}")

    script:
    """
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_7 \\
                                  --maf 0.05 \\
                                  --make-bed \\
                                  --out HapMap_3_r3_8
    """
}

process HWE_CHECKING {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path maf_filering_output

    output:
    path("plink.hwe")

    script:
    """
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_8 --hardy
    """
}

process HWE_CHECKING_PLOT {
    container "rocker/r-base:4.5.2"

    input:
    path hwe_check_output
    path hwe_plot_script

    output:
    path("histhwe.pdf"),   emit: normal_plot
    path("histhwe_below_theshold.pdf"), emit: zoom_plot

    script:
    """
    awk '{ if (\$9 <0.00001) print \$0 }' plink.hwe>plinkzoomhwe.hwe
    Rscript --no-save $hwe_plot_script
    """
}

process HWE_HANDLE_CONTROL {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path maf_filering_output

    output:
    path("HapMap_hwe_filter_step1.{bed,bim,fam}")

    script:
    """
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_8 \\
                                  --hwe 1e-6 \\
                                  --make-bed \\
                                  --out HapMap_hwe_filter_step1
    """ 
}

process HWE_HANDLE_CASES {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path first_hwe_filter

    output:
    path("HapMap_3_r3_9.{bed,bim,fam}")

    script:
    """
    /usr/lib/debian-med/bin/plink --bfile HapMap_hwe_filter_step1 \\
                                  --hwe 1e-10 \\
                                  --hwe-all \\
                                  --make-bed \\
                                  --out HapMap_3_r3_9
    """ 
}

process HETEROZYGOSITY_CHECK {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path hwe_filter_output
    path inversion_file

    output:
    path("indepSNP.prune.in"), emit: prune
    path("R_check.het"), emit: het

    script:
    """
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_9 \\
                                  --exclude inversion.txt \\
                                  --range \\
                                  --indep-pairwise 50 5 0.2 \\
                                  --out indepSNP

    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_9 \\
                                  --extract indepSNP.prune.in \\
                                  --het \\
                                  --out R_check
    """
}

process HETEROZYGOSITY_PLOT_DISTRIBUTION {
    container "rocker/r-base:4.5.2"

    input:
    path het_check
    path het_plot_rate_script
    path het_outlier_script

    output:
    path("heterozygosity.pdf"), emit: rate 
    path("fail-het-qc.txt"), emit: outlier

    script:
    """
    Rscript --no-save $het_plot_rate_script
    Rscript --no-save $het_outlier_script
    """
}

process HETEROZYGOSITY_DROP_FAILED_SAMPLES {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path hwe_filter_output
    path het_outlier_list 

    output:
    path("HapMap_3_r3_10.{bed,bim,fam}")

    script:
    """
    sed 's/"// g' $het_outlier_list| awk '{print\$1, \$2}'> het_fail_ind.txt
    
    /usr/lib/debian-med/bin/plink --bfile HapMap_3_r3_9 \\
                                  --remove het_fail_ind.txt \\
                                  --make-bed \\
                                  --out HapMap_3_r3_10
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

    // Step 1: Filtering missingness
    FILTERING_MISSINGNESS(input_files_ch)

    // Step 2: Check sex discrepancy
    CHECK_SEX_DISCREPANCY(FILTERING_MISSINGNESS.out.second_filter_individual)

    visualize_sex_script = channel.fromPath("${projectDir}/1_QC_GWAS/gender_check.R")
    SEX_DISCREPANCY_VISUALIZE(
        CHECK_SEX_DISCREPANCY.out,
        visualize_sex_script
    )

    SEX_DISCREPANCY_HANDLE_BY_REMOVE(
        CHECK_SEX_DISCREPANCY.out,
        FILTERING_MISSINGNESS.out.second_filter_individual
    )

    SEX_DISCREPANCY_HANDLE_BY_IMPUTE(
        FILTERING_MISSINGNESS.out.second_filter_individual
    )

    // Step 3: Filtering by MAF
    // Generate a bfile with autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF).
    MAF_AUTOSOMAL_SELECTION(SEX_DISCREPANCY_HANDLE_BY_IMPUTE.out)

    maf_check_script = channel.fromPath("${projectDir}/1_QC_GWAS/MAF_check.R")
    MAF_PLOT_DISTRIBUTION(
        MAF_AUTOSOMAL_SELECTION.out.check_frq, 
        maf_check_script
    )
    
    MAF_FILTERING(MAF_AUTOSOMAL_SELECTION.out.plink_files)

    // Step 4: Filtering by Hardy-Weinberg equilibrium (HWE)
    HWE_CHECKING(MAF_FILTERING.out)

    hwe_plot_script = channel.fromPath("${projectDir}/1_QC_GWAS/hwe.R")
    HWE_CHECKING_PLOT(
        HWE_CHECKING.out,
        hwe_plot_script 
    )

    HWE_HANDLE_CONTROL(MAF_FILTERING.out)
    HWE_HANDLE_CASES(HWE_HANDLE_CONTROL.out)
    
    // Step 5: Heterozygosity
    // Excludes individuals with high or low heterozygosity rates
    inversion_file = channel.fromPath("./1_QC_GWAS/inversion.txt")
    HETEROZYGOSITY_CHECK(HWE_HANDLE_CASES.out, inversion_file)

    het_plot_script = channel.fromPath("./1_QC_GWAS/check_heterozygosity_rate.R")
    het_outlier_script = channel.fromPath("./1_QC_GWAS/heterozygosity_outliers_list.R")
    HETEROZYGOSITY_PLOT_DISTRIBUTION(
        HETEROZYGOSITY_CHECK.out.het,
        het_plot_script,
        het_outlier_script
    )

    HETEROZYGOSITY_DROP_FAILED_SAMPLES(
        HWE_HANDLE_CASES.out,
        HETEROZYGOSITY_PLOT_DISTRIBUTION.out.outlier
    )
}