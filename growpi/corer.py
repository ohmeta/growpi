#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import textwrap
import shutil

import growpi

WORKFLOWS = ["grid_all", "all"]


def run_snakemake(args, unknown, snakefile, workflow):
    conf = growpi.parse_yaml(args.config)

    if not os.path.exists(conf["params"]["samples"]):
        print("Please specific samples list on init step or change config.yaml manualy")
        sys.exit(1)

    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--configfile",
        args.config,
        "--cores",
        str(args.cores),
        "--until",
        args.task
    ] + unknown

    if "--touch" in unknown:
        pass
    elif args.conda_create_envs_only:
        cmd += ["--use-conda", "--conda-create-envs-only"]
        if args.conda_prefix is not None:
            cmd += ["--conda-prefix", args.conda_prefix]
    else:
        cmd += [
            "--keep-going",
            "--printshellcmds",
            "--reason",
        ]

        if args.use_conda:
            cmd += ["--use-conda"]
            if args.conda_prefix is not None:
                cmd += ["--conda-prefix", args.conda_prefix]

        if args.list:
            cmd += ["--list"]
        elif args.run_local:
            cmd += ["--local-cores", str(args.local_cores),
                    "--jobs", str(args.jobs)]
        elif args.run_remote:
            cmd += ["--profile", args.profile,
                    "--local-cores", str(args.local_cores),
                    "--jobs", str(args.jobs)]
        elif args.debug:
            cmd += ["--debug-dag"]
        else:
            cmd += ["--dry-run"]

        if args.dry_run and ("--dry-run" not in cmd):
            cmd += ["--dry-run"]

    cmd_str = " ".join(cmd).strip()
    print("Running growpi %s:\n%s" % (workflow, cmd_str))

    env = os.environ.copy()
    proc = subprocess.Popen(
        cmd_str,
        shell=True,
        stdout=sys.stdout,
        stderr=sys.stderr,
        env=env,
    )
    proc.communicate()

    print(f'''\nReal running cmd:\n{cmd_str}''')


def update_env_path(conf, work_dir):
    for env_name in conf["envs"]:
        conf["envs"][env_name] = os.path.join(os.path.abspath(work_dir), conf["envs"][env_name])
    return conf


def update_config_tools(conf, begin, trimmer, rmhoster, profilers):
    conf["params"]["simulate"]["do"] = False
    conf["params"]["begin"] = begin

    for trimmer_ in ["sickle", "fastp"]:
        if trimmer_ == trimmer:
            conf["params"]["trimming"][trimmer_]["do"] = True
        else:
            conf["params"]["trimming"][trimmer_]["do"] = False

    for rmhoster_ in ["bwa", "bowtie2", "minimap2", "kraken2", "kneaddata"]:
        if rmhoster_ == rmhoster:
            conf["params"]["rmhost"][rmhoster_]["do"] = True
        else:
            conf["params"]["rmhost"][rmhoster_]["do"] = False

    for profiler_ in profilers:
        if profiler_ in ["irep", "grid", "smeg", "demic"]:
            if profiler_ in conf["params"]["profiling"]:
                conf["params"]["profiling"][profiler_]["do"] = True
            else:
                print(f"Currently growpi don't support {profiler_}")
        else:
            if profiler_ in conf["params"]["profiling"]:
                conf["params"]["profiling"][profiler_]["do"] = False
            else:
                print(f"Currently growpi don't support {profiler_}")

    if begin == "simulate":
        conf["params"]["simulate"]["do"] = True
    elif begin == "rmhost":
        conf["params"]["trimming"][trimmer]["do"] = False
    elif begin == "profiling":
        conf["params"]["raw"]["save_reads"] = True
        conf["params"]["raw"]["fastqc"]["do"] = False
        conf["params"]["qcreport"]["do"] = False

        conf["params"]["trimming"][trimmer]["do"] = False
        conf["params"]["rmhost"][rmhoster]["do"] = False
    return conf


def init(args, unknown):
    if args.workdir:
        project = growpi.metaconfig(args.workdir)
        print(project.__str__())
        project.create_dirs()
        conf = project.get_config()

        conf = update_config_tools(
            conf, args.begin, args.trimmer, args.rmhoster, args.profiler
        )
        conf = update_env_path(conf, args.workdir)

        if args.samples:
            conf["params"]["samples"] = args.samples
        else:
            print("Please supply samples table")
            sys.exit(-1)

        growpi.update_config(
            project.config_file, project.new_config_file, conf, remove=False
        )
    else:
        print("Please supply a workdir!")
        sys.exit(-1)


def run_wf(args, unknown):
    snakefile = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                             "Snakefile")
    run_snakemake(args, unknown, snakefile)


def main():
    banner = """

    ▄████  ██▀███   ▒█████   █     █░ ██▓███   ██▓
    ██▒ ▀█▒▓██ ▒ ██▒▒██▒  ██▒▓█░ █ ░█░▓██░  ██▒▓██▒
    ▒██░▄▄▄░▓██ ░▄█ ▒▒██░  ██▒▒█░ █ ░█ ▓██░ ██▓▒▒██▒
    ░▓█  ██▓▒██▀▀█▄  ▒██   ██░░█░ █ ░█ ▒██▄█▓▒ ▒░██░
    ░▒▓███▀▒░██▓ ▒██▒░ ████▓▒░░░██▒██▓ ▒██▒ ░  ░░██░
    ░▒   ▒ ░ ▒▓ ░▒▓░░ ▒░▒░▒░ ░ ▓░▒ ▒  ▒▓▒░ ░  ░░▓  
    ░   ░   ░▒ ░ ▒░  ░ ▒ ▒░   ▒ ░ ░  ░▒ ░      ▒ ░
    ░ ░   ░   ░░   ░ ░ ░ ░ ▒    ░   ░  ░░        ▒ ░
        ░    ░         ░ ░      ░              ░  
 
           Omics for All, Open Source for All

    Microbial community grow rate profiling pipeline

"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(banner),
        prog="growpi",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print software version and exit",
    )

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "-d",
        "--workdir",
        metavar="WORKDIR",
        type=str,
        default="./",
        help="project workdir",
    )
    common_parser.add_argument(
        "--check-samples",
        dest="check_samples",
        default=False,
        action="store_true",
        help="check samples, default: False",
    )

    run_parser = argparse.ArgumentParser(add_help=False)
    run_parser.add_argument(
        "--config",
        type=str,
        default="./config.yaml",
        help="config.yaml",
    )
    run_parser.add_argument(
        "--profile",
        type=str,
        default="./profiles/slurm",
        help="cluster profile name",
    )
    run_parser.add_argument(
        "--cores",
        type=int,
        default=240,
        help="all job cores, available on '--run-local'")
    run_parser.add_argument(
        "--local-cores",
        type=int,
        dest="local_cores",
        default=8,
        help="local job cores, available on '--run-remote'")
    run_parser.add_argument(
        "--jobs",
        type=int,
        default=30,
        help="cluster job numbers, available on '--run-remote'")
    run_parser.add_argument(
        "--list",
        default=False,
        action="store_true",
        help="list pipeline rules",
    )
    run_parser.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="debug pipeline",
    )
    run_parser.add_argument(
        "--dry-run",
        default=False,
        dest="dry_run",
        action="store_true",
        help="dry run pipeline",
    )
    run_parser.add_argument(
        "--run-local",
        default=False,
        action="store_true",
        help="run pipeline on local computer",
    )
    run_parser.add_argument(
        "--run-remote",
        default=False,
        action="store_true",
        help="run pipeline on remote cluster",
    )
    run_parser.add_argument(
        "--cluster-engine",
        default="slurm",
        choices=["slurm", "sge", "lsf", "pbs-torque"],
        help="cluster workflow manager engine, support slurm(sbatch) and sge(qsub)"
    )
    run_parser.add_argument("--wait", type=int, default=60, help="wait given seconds")
    run_parser.add_argument(
        "--use-conda",
        default=False,
        dest="use_conda",
        action="store_true",
        help="use conda environment",
    )
    run_parser.add_argument(
        "--conda-prefix",
        default="~/.conda/envs",
        dest="conda_prefix",
        help="conda environment prefix",
    )
    run_parser.add_argument(
        "--conda-create-envs-only",
        default=False,
        dest="conda_create_envs_only",
        action="store_true",
        help="conda create environments only",
    )

    subparsers = parser.add_subparsers(title="available subcommands", metavar="")
    parser_init = subparsers.add_parser(
        "init",
        formatter_class=growpi.custom_help_formatter,
        parents=[common_parser],
        prog="growpi init",
        help="init project",
    )
    parser_run_wf = subparsers.add_parser(
        "run_wf",
        formatter_class=growpi.custom_help_formatter,
        parents=[common_parser, run_parser],
        prog="growpi run_wf",
        help="metagenome grow rate profiling pipeline",
    )

    parser_init.add_argument(
        "-s",
        "--samples",
        type=str,
        default=None,
        help="""desired input:
samples list, tsv format required.

if begin from trimming, rmhost, or assembly:
    if it is fastq:
        the header is: [sample_id, fq1, fq2]
    if it is sra:
        the header is: [sample_id, sra]

if begin from simulate:
        the header is: [id, genome, abundance, reads_num, model]

""",
    )
    parser_init.add_argument(
        "-b",
        "--begin",
        type=str,
        default="trimming",
        choices=["simulate", "trimming", "rmhost", "profiling"],
        help="pipeline starting point",
    )
    parser_init.add_argument(
        "--trimmer",
        type=str,
        default="fastp",
        required=False,
        choices=["oas1", "sickle", "fastp"],
        help="which trimmer used",
    )
    parser_init.add_argument(
        "--rmhoster",
        type=str,
        default="bowtie2",
        required=False,
        choices=["bwa", "bowtie2", "minimap2", "kraken2", "kneaddata"],
        help="which rmhoster used",
    )
    parser_init.add_argument(
        "--profiler",
        nargs="+",
        default=["grid"],
        required=False,
        choices=["irep", "grid", "smeg", "demic"],
        help="which profiler used",
    )
 
    parser_init.set_defaults(func=init)

    parser_run_wf.add_argument(
        "task",
        metavar="TASK",
        nargs="?",
        type=str,
        default="all",
        choices=WORKFLOWS,
        help="pipeline end point. Allowed values are " + ", ".join(WORKFLOWS),
    )
    parser_run_wf.set_defaults(func=run_wf)

    args, unknown = parser.parse_known_args()

    try:
        if args.version:
            print("growpi version %s" % growpi.__version__)
            sys.exit(0)
        args.func(args, unknown)
    except AttributeError as e:
        print(e)
        parser.print_help()


if __name__ == "__main__":
    main()
