configfile: "dsprint/config.json"

rule csq:
    threads: 4
    input:
        config['csq']['input_file']
    output:
        output=expand(
            "{output_prefix}{output}{output_suffix}",
            output_prefix=config['csq']['output_prefix'],
            output=list(range(1, 23)) + ['X', 'Y'],
            output_suffix=config['csq']['output_suffix']
        )
    script:
        "scripts/1.parse_ExAC/ExAC_parser.py"
