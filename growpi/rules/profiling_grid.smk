def profiling_input_with_short_reads(wildcards, have_single=False):
    if RMHOST_DO:
        return get_reads(wildcards, "rmhost", have_single)
    if TRIMMING_DO:
        return get_reads(wildcards, "trimming", have_single)
    else:
        return get_reads(wildcards, "raw", have_single)


rule profiling_prepare_reads:
    input:
        profiling_input_with_short_reads
    output:
        temp(os.path.join(config["output"]["profiling"], "fastq/{sample}/{sample}.fastq.gz"))
    shell:
        '''
        seqtk mergepe {input} | \
        pigz -fc > {output}
        '''


if config["params"]["profiling"]["grid"]["do_single"]:
    rule profiling_grid_single:
        input:
            os.path.join(config["output"]["profiling"], "fastq/grid/{sample}.fastq.gz")
        output:
            txt = os.path.join(config["output"]["profiling"], "profile/grid_single/{sample}/{sample}.GRiD.txt")
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/grid_single/{sample}.grid.benchmark.txt")
        log:
            os.path.join(config["output"]["profiling"], "logs/grid_single/{sample}.grid.log")
        conda:
            config["envs"]["grid"]
        params:
            sample_id = "{sample}",
            database = config["params"]["profiling"]["grid"]["database"],
            coverage_cutoff = config["params"]["profiling"]["grid"]["coverage_cutoff"],
            reassignment = "-p" if config["params"]["profiling"]["grid"]["reassignment"] else "",
            out_dir = os.path.join(config["output"]["profiling"], "profile/grid_single/{sample}"),
            fq_dir = os.path.join(config["output"]["profiling"], "fastq/{sample}")
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            grid single \
            -r {params.fq_dir} \
            -e fastq.gz \
            -o {params.out_dir} \
            {params.reassignment} \
            -c {params.coverage_cutoff} \
            -d {params.database} \
            -n {threads} \
            >{log} 2>&1
            '''


    rule profiling_grid_single_merge:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/grid_single/{sample}/{sample}.GRiD.txt"),
                sample=SAMPLES.index.unique())
        output:
            grid = os.path.join(config["output"]["profiling"], "report/grid_single/grid.merged.tsv"),
            grid_reliable = os.path.join(config["output"]["profiling"], "report/grid_single/grid.merged.reliable.tsv")
        threads:
            config["params"]["profiling"]["threads"]
        run:
            import os
            import pandas as pd

            def parse_grid(txt_file):
                try:
                    if os.path.exists(txt_file):
                        df = pd.read_csv(txt_file, sep="\t")
                        if not df.empty:
                            sample_id = os.path.basename(txt_file).split(".")[0]
                            df_reliable = df.query('Species_heterogeneity < 0.3')
                            key = ["Genome", "GRiD"]
                            df = df.loc[:, key].rename(columns={"GRiD": sample_id}).set_index("Genome")
                            df_reliable = df.loc[:, key].rename(columns={"GRiD": sample_id}).set_index("Genome")
                            return df, df_reliable
                        else:
                            return None, None
                    else:
                        print("%s is not exists" % txt_file)
                        return None, None
                except pd.io.common.EmptyDataError:
                    print("%s is empty, please check" % txt_file)
                    return None, None

            growpi.merge(input, threads, parse_grid,
                         output_1=output.grid, output_2=output.grid_reliable)


    rule profiling_grid_single_all:
        input:
            os.path.join(config["output"]["profiling"], "report/grid_single/grid.merged.tsv"),
            os.path.join(config["output"]["profiling"], "report/grid_single/grid.merged.reliable.tsv")


else:
    rule profiling_grid_single_all:
        input:


if config["params"]["profiling"]["grid"]["do_multi"]:
    rule profiling_grid_multi:
        input:
            os.path.join(config["output"]["profiling"], "fastq/grid/{sample}.fastq.gz")
        output:
            txt = os.path.join(config["output"]["profiling"], "profile/grid_multi/{sample}/{sample}.GRiD.txt")
        benchmark:
            os.path.join(config["output"]["profiling"], "benchmark/grid_multi/{sample}.grid.benchmark.txt")
        log:
            os.path.join(config["output"]["profiling"], "logs/grid_multi/{sample}.grid.log")
        conda:
            config["envs"]["grid"]
        params:
            sample_id = "{sample}",
            database = config["params"]["profiling"]["grid"]["database"],
            coverage_cutoff = config["params"]["profiling"]["grid"]["coverage_cutoff"],
            reassignment = "-p" if config["params"]["profiling"]["grid"]["reassignment"] else "",
            out_dir = os.path.join(config["output"]["profiling"], "profile/grid_multi/{sample}"),
            fq_dir = os.path.join(config["output"]["profiling"], "fastq/{sample}")
        threads:
            config["params"]["profiling"]["threads"]
        shell:
            '''
            grid multiplex \
            -r {params.fq_dir} \
            -e fastq.gz \
            -o {params.out_dir} \
            {params.reassignment} \
            -c {params.coverage_cutoff} \
            -d {params.database} \
            -n {threads} \
            >{log} 2>&1
            '''


    rule profiling_grid_multi_merge:
        input:
            expand(os.path.join(
                config["output"]["profiling"],
                "profile/grid_multi/{sample}/{sample}.GRiD.txt"),
                sample=SAMPLES.index.unique())
        output:
            grid = os.path.join(config["output"]["profiling"], "report/grid_multi/grid.merged.tsv"),
            grid_reliable = os.path.join(config["output"]["profiling"], "report/grid_multi/grid.merged.reliable.tsv")
        threads:
            config["params"]["profiling"]["threads"]
        run:
            import os
            import pandas as pd

            def parse_grid(txt_file):
                try:
                    if os.path.exists(txt_file):
                        df = pd.read_csv(txt_file, sep="\t")
                        if not df.empty:
                            sample_id = os.path.basename(txt_file).split(".")[0]
                            df_reliable = df.query('Species_heterogeneity < 0.3')
                            key = ["Genome", "GRiD"]
                            df = df.loc[:, key].rename(columns={"GRiD": sample_id}).set_index("Genome")
                            df_reliable = df.loc[:, key].rename(columns={"GRiD": sample_id}).set_index("Genome")
                            return df, df_reliable
                        else:
                            return None, None
                    else:
                        print("%s is not exists" % txt_file)
                        return None, None
                except pd.io.common.EmptyDataError:
                    print("%s is empty, please check" % txt_file)
                    return None, None

            growpi.merge(input, threads, parse_grid,
                         output_1=output.grid, output_2=output.grid_reliable)


    rule profiling_grid_multi_all:
        input:
            os.path.join(config["output"]["profiling"], "report/grid_multi/grid.merged.tsv"),
            os.path.join(config["output"]["profiling"], "report/grid_multi/grid.merged.reliable.tsv")


else:
    rule profiling_grid_multi_all:
        input: