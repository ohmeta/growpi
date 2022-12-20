rule grid:
    input:
        fq1 = lambda wildcards: metapi.get_reads(SAMPLES, wildcards, "fq1")[0],
        fq2 = lambda wildcards: metapi.get_reads(SAMPLES, wildcards, "fq2")[0]
    output:
        fq_dir = directory(temp(os.path.join(config["output"]["grid"],
                                             "{sample}.fq.temp"))),
        txt = os.path.join(config["output"]["grid"],
                           "{sample}.grid.out/{sample}.GRiD.txt")
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


rule grid_merge:
    input:
        expand(os.path.join(config["output"]["grid"],
                            "{sample}.grid.out/{sample}.GRiD.txt"),
               sample=SAMPLES.index.unique())
    output:
        grid = os.path.join(config["output"]["grid"], "grid.merged.tsv"),
        grid_reliable = os.path.join(config["output"]["grid"], "grid.merged.reliable.tsv")
    threads:
        config["params"]["grid"]["threads"]
    run:
        import os
        import pandas as pd
        import concurrent.futures

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


        def merge(files, workers, func, **kwargs):
            df_1_list = []
            df_2_list = []

            with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                for df_1, df_2 in executor.map(func, files):
                    if df_1 is not None:
                        df_1_list.append(df_1)
                    if df_2 is not None:
                        df_2_list.append(df_2)

            df_1_ = pd.concat(df_1_list, axis=1).fillna(1)
            df_2_ = pd.concat(df_2_list, axis=1).fillna(1)

            if "output_1" in kwargs:
                df_1_.reset_index().to_csv(kwargs["output_1"], sep="\t", index=False)
            if "output_2" in kwargs:
                df_2_.reset_index().to_csv(kwargs["output_2"], sep="\t", index=False)

            return df_1_, df_2_


        merge(input, threads, parse_grid, output_1=output.grid, output_2=output.grid_reliable)


rule grid_all:
    input:
        os.path.join(config["output"]["grid"], "grid.merged.tsv"),
        os.path.join(config["output"]["grid"], "grid.merged.reliable.tsv")
