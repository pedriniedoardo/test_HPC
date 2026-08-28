rule testRunPixi:
    '''
    This is a test rule for pixi
    '''
    output:
        out_file = config["dir_output"] + config["test_out"]
    conda: 
        config["pixi_env"]
    log:
        'logs/testRunPixi.log'
    benchmark:
        'benchmarks/testRunPixi.log'
    shell:
        '''
        pixi list > {output.out_file}
        '''

rule testRunR:
    '''
    This is a test rule for R and renv integration
    '''
    output:
        out_file = config["dir_output"] + config["test_outR"]
    conda: 
        config["pixi_env"]
    log:
        'logs/testRunR.log'
    benchmark:
        'benchmarks/testRunR.log'
    script:
        "../../scr/snakemake/000_testR_snakemake.R"
