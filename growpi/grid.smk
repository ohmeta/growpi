rule grid:
    input:
        fq1 = lambda wildcards: metapi.get_reads(SAMPLES, wildcards, "fq1")[0],
        fq2 = lambda wildcards: metapi.get_reads(SAMPLES, wildcards, "fq2")[0]
    output:
        fq_dir = directory(temp(os.path.join(config["output"]["grid"],
                                             "{sample}.fq.temp"))),
        done = os.path.join(config["output"]["grid"], "{sample}.grid.out/done")
    log:
        os.path.join(config["output"]["grid"], "logs/{sample}.grid.log")
    conda:
        config["envs"]["bioenv2"]
    params:
        id = "{sample}",
        database = config["params"]["grid"]["database"],
        coverage_cutoff = config["params"]["grid"]["coverage_cutoff"],
        reassignment = "-p" if config["params"]["grid"]["reassignment"] \
            else "",
        out_dir = os.path.join(config["output"]["grid"], "{sample}.grid.out")

    threads:
        config["params"]["grid"]["threads"]
    shell:
        '''
        mkdir -p {output.fq_dir}
        seqtk mergepe {input.fq1} {input.fq2} | \
        pigz -c > {output.fq_dir}/{params.id}.fastq.gz

        grid multiplex \
        -r {output.fq_dir} \
        -e fastq.gz \
        -o {params.out_dir} \
        {params.reassignment} \
        -c {params.coverage_cutoff} \
        -d {params.database} \
        -n {threads} > {log} 2>&1

        touch {output.done}
        '''


rule grid_all:
    input:
        expand(os.path.join(config["output"]["grid"],
                            "{sample}.grid.out/done"),
               sample=SAMPLES.index.unique())
