rule csq:
    input:
        "/media/vineetb/t5-vineetb/dsprint/ExAC.r0.3.sites.vep.vcf"
    output:
        expand("/media/vineetb/t5-vineetb/dsprint/out/parsed_chrom_{c}.csv", c=[1,2,4])
    script:
        "scripts/1.parse_ExAC/ExAC_parser.py"
