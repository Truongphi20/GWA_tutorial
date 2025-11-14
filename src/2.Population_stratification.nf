process OKGP_DOWNLOAD {

    scratch "${params.store_dir}"

    output:
    path("ALL.2of4intersection.20100804.genotypes.vcf.gz")

    script:
    """
    wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz
    """
}


process OKGP_CONVERT_BFILE {
    container "biocontainers/plink:v1.07dfsg-2-deb_cv1"

    input:
    path okgp_vcf

    output:
    path("ALL.2of4intersection.20100804.genotypes.{bed,bim,fam}")

    script:
    """
    wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz

    /usr/lib/debian-med/bin/plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz \\
                                  --make-bed \\
                                  --out ALL.2of4intersection.20100804.genotypes
    """
}


workflow POP_STRATIFICATION {
    // Step 1: Prepare 1KGP data 
    OKGP_DOWNLOAD()
    OKGP_CONVERT_BFILE(OKGP_DOWNLOAD.out)
}