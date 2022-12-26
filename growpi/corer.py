#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import textwrap
import shutil
import metapi


WORKFLOWS = ["grid_all", "all"]


def run_snakemake(args, unknown, snakefile):
    conf = metapi.parse_yaml(args.config)

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
    ] + unknown

    if args.conda_create_envs_only:
        cmd += ["--use-conda", "--conda-create-envs-only"]
    else:
        cmd += [
            "--rerun-incomplete",
            "--keep-going",
            "--printshellcmds",
            "--reason",
            "--until",
            args.task,
        ]

        if args.use_conda:
            cmd += ["--use-conda"]

        if args.list:
            cmd += ["--list"]
        elif args.run:
            cmd += [""]
        elif args.debug:
            cmd += ["--debug-dag", "--dry-run"]
        elif args.dry_run:
            cmd += ["--dry-run"]
        elif args.qsub:
            cmd += [
                "--cluster-config",
                args.cluster,
                "--jobs",
                str(args.jobs),
                "--latency-wait",
                str(args.wait),
                '--cluster "qsub -S /bin/bash -cwd \
                -q {cluster.queue} -P {cluster.project} \
                -l vf={cluster.mem},p={cluster.cores} \
                -binding linear:{cluster.cores} \
                -o {cluster.output} -e {cluster.error}"',
            ]

    cmd_str = " ".join(cmd).strip()
    print(f"Running growpi:\n{cmd_str}")

    env = os.environ.copy()
    proc = subprocess.Popen(
        cmd_str, shell=True, stdout=sys.stdout, stderr=sys.stderr, env=env,
    )
    proc.communicate()


class growpi_config:
    sub_dirs = ["envs", "results", "logs/grid", "logs/grid_merge"]

    def __init__(self, work_dir):
        self.work_dir = os.path.realpath(work_dir)
        self.growpi_dir = os.path.dirname(os.path.abspath(__file__))

        self.config_file = os.path.join(self.growpi_dir, "config", "config.yaml")
        self.cluster_file = os.path.join(self.growpi_dir, "config", "cluster.yaml")
        self.snake_file = os.path.join(self.growpi_dir, "Snakefile")
        self.envs_dir = os.path.join(self.growpi_dir, "envs")

        self.new_config_file = os.path.join(self.work_dir, "config.yaml")
        self.new_cluster_file = os.path.join(self.work_dir, "cluster.yaml")

    def __str__(self):
        message = """
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

    Microbial commnuity grow rate profiling pipeline

Thanks for using growpi.

A metagenomics project has been created at %s


if you want to create fresh conda environments:

        growpi run_wf --conda_create_envs_only

if you have environments:

        growpi --help
        growpi init --help
        growpi run_wf --help
""" % (
            self.work_dir
        )

        return message

    def create_dirs(self):
        """
        create project directory
        """
        if not os.path.exists(self.work_dir):
            os.mkdir(self.work_dir)

        for sub_dir in growpi_config.sub_dirs:
            os.makedirs(os.path.join(self.work_dir, sub_dir), exist_ok=True)

        for i in os.listdir(self.envs_dir):
            shutil.copyfile(
                os.path.join(self.envs_dir, i),
                os.path.join(self.work_dir, os.path.join("envs", i)),
            )

    def get_config(self):
        """
        get default configuration
        """
        config = metapi.parse_yaml(self.config_file)
        cluster = metapi.parse_yaml(self.cluster_file)
        return (config, cluster)


def init(args, unknown):
    if args.workdir:
        project = growpi_config(args.workdir)
        print(project.__str__())
        project.create_dirs()

        conf, cluster = project.get_config()

        conf["envs"]["bioenv2"] = os.path.join(
            os.path.realpath(args.workdir), "envs/bioenv2.yaml"
        )

        if args.samples:
            conf["params"]["samples"] = args.samples

        metapi.update_config(
            project.config_file, project.new_config_file, conf, remove=False
        )
        metapi.update_config(
            project.cluster_file, project.new_cluster_file, cluster, remove=False
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

      ██████  ██████   ██████  ██     ██ ██████  ██ 
     ██       ██   ██ ██    ██ ██     ██ ██   ██ ██ 
     ██   ███ ██████  ██    ██ ██  █  ██ ██████  ██ 
     ██    ██ ██   ██ ██    ██ ██ ███ ██ ██      ██ 
      ██████  ██   ██  ██████   ███ ███  ██      ██ 

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
        help="project workdir, default: ./",
    )

    run_parser = argparse.ArgumentParser(add_help=False)
    run_parser.add_argument(
        "--config",
        type=str,
        default="./config.yaml",
        help="config.yaml, default: ./config.yaml",
    )
    run_parser.add_argument(
        "--cluster",
        type=str,
        default="./cluster.yaml",
        help="cluster.yaml, default: ./cluster.yaml",
    )
    run_parser.add_argument(
        "--cores", type=int, default=8, help="CPU cores, default: 8"
    )
    run_parser.add_argument(
        "--jobs", type=int, default=80, help="qsub job numbers, default: 80"
    )
    run_parser.add_argument(
        "--list", default=False, action="store_true", help="list pipeline rules",
    )
    run_parser.add_argument(
        "--run", default=False, action="store_true", help="run pipeline",
    )
    run_parser.add_argument(
        "--debug", default=False, action="store_true", help="debug pipeline",
    )
    run_parser.add_argument(
        "--dry_run", default=False, action="store_true", help="dry run pipeline",
    )
    run_parser.add_argument(
        "--qsub", default=False, action="store_true", help="qsub pipeline",
    )
    run_parser.add_argument(
        "--wait", type=int, default=60, help="wait given seconds, default: 60"
    )
    run_parser.add_argument(
        "--use_conda", default=False, action="store_true", help="use conda environment"
    )
    run_parser.add_argument(
        "--conda_create_envs_only",
        default=False,
        action="store_true",
        help="conda create environments only",
    )

    subparsers = parser.add_subparsers(title="available subcommands", metavar="")
    parser_init = subparsers.add_parser(
        "init", parents=[common_parser], prog="growpi init", help="init project",
    )
    parser_run_wf = subparsers.add_parser(
        "run_wf",
        parents=[common_parser, run_parser],
        prog="growpi run_wf",
        help="computing grow rate",
    )

    parser_init.add_argument(
        "-s",
        "--samples",
        type=str,
        default=None,
        help="""desired input:
samples list, tsv format required.
        the header is: [id, fq1, fq2]
""",
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
