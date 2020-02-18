configfile: "dsprint/config.json"
workdir: "/home/vineetb/git_checkouts/dsprint"
threads: 4

rule all:
    input:
        expand(
            "{output_prefix}{chromosome}{output_suffix}",
            output_prefix=config['csq']['output_prefix'],
            chromosome=list(range(1,23) )+ ['X','Y'],
            output_suffix=config['csq']['output_suffix']
        )

# Jobs run much faster if they don't have to independently gunzip the file
rule gunzip:
    input:
        config['csq']['input_file']
    output:
        config['csq']['input_file'][:-3]
    shell:
        "gunzip --keep {input}"

rule csq:
    input:
        config['csq']['input_file'][:-3]
    output:
        config['csq']['output_prefix'] + "{chromosome}" + config['csq']['output_suffix']
    wildcard_constraints:
        chromosome="\d+|X|Y"
    script:
        "scripts/1.parse_ExAC/ExAC_parser.py"
