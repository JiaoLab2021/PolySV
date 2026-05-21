#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"


import os
import argparse
import sys
import logging
import shutil


# environment variables
env_path = {'PATH': os.environ.get('PATH')}
# include
base_dir = os.path.dirname(os.path.abspath(__file__))
target_include = os.path.join(base_dir, "include")
if os.path.exists(target_include):
    sys.path.insert(0, target_include)
    for root, dirs, files in os.walk(target_include):
        sys.path.insert(0, root)

# from saveResult import MySaveResult
from make_dir import MyMkdir
from align import MyAlign
from call import MyCall
from genome import MyGenome
from vcf2convert import MyVcf2Convert
from merge import MyMerge
from refine import MyRefine
from vcf_record import MyVcfRecord


# get parser
class MyParser:
    def __init__(self):
        # logging
        self.set_logging()
        # parameter
        self.set_parameter()
        self.args = self.parser.parse_args()
        # path
        self.parameter_path()

    # logging
    def set_logging(self):
        # Configure the logging
        self.logger = logging.getLogger('polysv')  # create log object
        self.logger.error(f"data: {__data__}")
        self.logger.error(f"version: {__version__}")
        self.logger.error(f"author: {__author__}")
        self.logger.error(f"If you encounter any issues related to the code, "
                          f"please don't hesitate to contact us via email at {__email__}.\n")
        formatter = logging.Formatter('[%(asctime)s] %(message)s')
        handler = logging.StreamHandler()  # output to the console
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

    # parameter
    def set_parameter(self):
        # parameters
        self.parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        # Input parameters
        self.required_input = self.parser.add_argument_group("Input arguments")

        self.required_input.add_argument(
            "-r", "--reference",
            dest="reference",
            help="reference genome file",
            default="refgenome.fa",
            required=True
        )
        self.required_input.add_argument(
            "-tgs", "--fq",
            dest="fq",
            help="third-generation technology, Suffix is fastq(fq).gz or fastq(fq)",
            default="sample.fastq.gz",
        )
        self.required_input.add_argument(
            "-tgsdt", "--fqdt",
            dest="fqdt",
            help="third-generation sequencing data types: ccs or clr or ont",
            default="ccs",
        )
        self.required_input.add_argument(
            "-genome", "--fa",
            dest="fa",
            help="query genome data, accepted formats: .fa/.fasta (optionally gzipped)",
            default="sample.fa",
        )
        self.required_input.add_argument(
            "-n", "--name",
            dest="samplename",
            help="sample name",
            default="sample",
            required=True
        )

        # Output parameters
        self.required_output = self.parser.add_argument_group("Output arguments")

        self.required_output.add_argument(
            "-o", "--output",
            dest="output",
            help="output file name",
            type=str,
            default="sample.vcf.gz",
            required=True
        )

        # Optional parameters
        self.optional = self.parser.add_argument_group("Optional arguments")

        self.optional.add_argument(
            "-w", "--workdir",
            dest="workdir",
            help="work directory",
            type=str,
            default="work_dir",
        )
        self.optional.add_argument(
            "-j", "--jobs",
            dest="jobs",
            help="number of jobs",
            type=int,
            default=2
        )
        self.optional.add_argument(
            "-t", "--threads",
            dest="threads",
            help="number of threads",
            type=int,
            default=12
        )
        self.optional.add_argument(
            "--restart",
            dest="restart",
            action="store_true",
            help="attempt to restart existing workflow based on the state of output files [False]",
        )
        self.optional.add_argument(
            "--remove",
            dest="remove",
            action="store_true",
            help="delete the generated intermediate files [False]",
        )

        # Algorithm parameters
        self.algorithm = self.parser.add_argument_group("Algorithm arguments")

        self.algorithm.add_argument(
            "--depth_correction_l_threshold",
            dest="l_depth_ratio",
            help="threshold for correction depth based on TGS alignment",
            type=float,
            default=0.75
        )
        self.algorithm.add_argument(
            "--depth_correction_a_threshold",
            dest="a_depth_ratio",
            help="threshold for correction depth based on genome alignment",
            type=float,
            default=0.25
        )
        self.algorithm.add_argument(
            "--del_breakpoint_window",
            dest="del_window",
            help="the size of the sliding window in the DEL correction algorithm (bp)",
            type=int,
            default=25
        )
        self.algorithm.add_argument(
            "--ins_breakpoint_window",
            dest="ins_window",
            help="the size of the sliding window in the INS correction algorithm (bp)",
            type=int,
            default=3
        )

    # path
    def parameter_path(self):
        # input
        self.jobs = self.args.jobs
        self.threads = self.args.threads
        self.samplename = self.args.samplename
        self.restart = self.args.restart
        self.remove = self.args.remove
        # Check if the parameter exists
        if not self.args.fq and not self.args.fa:
            raise ValueError("You must provide at least one of -tgs (fq) or -genome (fa).")
        if bool(self.args.fq) != bool(self.args.fqdt):
            raise ValueError("If you provide the -tgs (fq) argument, you must also "
                             "provide the -tgsdt (fqdt): ccs, clr, ont")
        # Continue to assign
        self.fq_path = os.path.abspath(self.args.fq) if self.args.fq else None
        self.fa_path = os.path.abspath(self.args.fa) if self.args.fa else None
        self.reference_path = os.path.abspath(self.args.reference)
        self.output_path = os.path.abspath(self.args.output)
        self.workdir_path = os.path.abspath(self.args.workdir)
        self.fqdt = self.args.fqdt if self.args.fq else None

        # algorithm
        self.l_depth_ratio = self.args.l_depth_ratio
        self.a_depth_ratio = self.args.a_depth_ratio
        self.del_window = self.args.del_window
        self.ins_window= self.args.ins_window

        # Check if the file exists
        if not os.path.exists(self.args.reference):
            self.logger.error(f"Error: {self.args.reference} No such file or directory.")
            sys.exit(1)
        if self.args.fq and not os.path.exists(self.args.fq):
            self.logger.error(f"Error: {self.args.fq} No such file or directory.")
            sys.exit(1)
        if self.args.fa and not os.path.exists(self.args.fa):
            self.logger.error(f"Error: {self.args.fa} No such file or directory.")
            sys.exit(1)


def main():
    # Get the parsed parameters
    parser = MyParser()
    # Start operation
    parser.logger.error(f'Running ...')
    # Set working directory
    parser.logger.error(f' ')
    parser.logger.error(f'##### Create a working directory #####')
    mkdirClass = MyMkdir()
    mkdirClass.make(parser.workdir_path, parser.restart)
    os.chdir(parser.workdir_path)
    parser.logger.error(f"Set working directory: {parser.workdir_path}")
    # Align
    parser.logger.error(f' ')
    parser.logger.error(f'##### Align #####')
    alignClass = MyAlign(parser.reference_path, parser.fq_path, parser.fa_path, parser.fqdt, parser.samplename,
                         parser.jobs, parser.threads, parser.workdir_path, parser.restart)
    align_path_dict = alignClass.run_align()
    # Call
    parser.logger.error(f' ')
    parser.logger.error(f'##### Call SV #####')
    callClass = MyCall(align_path_dict, parser.reference_path, parser.fq_path, parser.fa_path, parser.fqdt,
                        parser.samplename, parser.jobs, parser.threads, parser.workdir_path, parser.restart)
    call_path_dict = callClass.run_call()
    # Build refgenome index
    parser.logger.error(f' ')
    parser.logger.error(f'##### Build refgenome index #####')
    referenceClass = MyGenome(parser.reference_path)
    sequenceDict = referenceClass.build_index()
    # Convert
    parser.logger.error(f' ')
    parser.logger.error(f'##### VCF format conversion #####')
    vcf2convertClass = MyVcf2Convert(call_path_dict, sequenceDict, parser.reference_path, parser.fq_path,
                                     parser.fa_path, parser.fqdt, parser.samplename, parser.jobs, parser.threads,
                                     parser.workdir_path, parser.restart)
    vcf2convert_path_dict = vcf2convertClass.run_vcf2convert()
    # Merge
    parser.logger.error(f' ')
    parser.logger.error(f'##### Merge #####')
    mergeClass = MyMerge(vcf2convert_path_dict, sequenceDict, parser.reference_path, parser.fq_path, parser.fa_path,
                               parser.samplename, parser.jobs, parser.threads, parser.workdir_path, parser.restart)
    merge_path = mergeClass.run_merge()
    # Build VCF index
    vcflist_path = list(vcf2convert_path_dict.values())
    vcflistClass = MyVcfRecord(vcflist_path)
    vcfDict = vcflistClass.build_dict()
    # Refine
    command_line_args  = sys.argv
    command_line_str = " ".join(command_line_args ) + '\n'
    parser.logger.error(f' ')
    parser.logger.error(f'##### Refine #####')
    refineClass = MyRefine(merge_path, align_path_dict, vcfDict, sequenceDict, parser.reference_path, parser.fq_path,
                           parser.fa_path, parser.output_path, parser.samplename, parser.jobs, parser.threads,
                           parser.workdir_path, parser.restart, parser.l_depth_ratio, parser.a_depth_ratio,
                           parser.del_window, parser.ins_window, command_line_str)
    refineClass.run_refine()
    # Remove
    if parser.remove:
        parser.logger.error(f'Delete all intermediate files ...')
        shutil.rmtree(parser.workdir_path)
    parser.logger.error(f'Done ...')


if __name__ == '__main__':
    main()
