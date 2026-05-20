#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"


import os
import run_cmd
import logging
import sys
import concurrent.futures


# environment variables
env_path = {'PATH': os.environ.get('PATH')}
# include
include_dir = os.path.dirname(os.path.realpath(__file__))
for root, dirs, files in os.walk(os.path.join(include_dir, "include")):
    sys.path.append(root)

from make_dir import MyMkdir
from vcf_io import MyVCF


class MyDepth():
    # Flag whether global log handler has been added
    logger_initialized = False
    def __init__(self, samplename, lbam_path, abam_path, depth_path, jobs, threads, restart, sequenceDict):
        # If the log handler has not been added yet, add
        if not MyDepth.logger_initialized:
            self.logger = logging.getLogger(__name__)
            formatter = logging.Formatter('[%(asctime)s] %(message)s')  # add date and time formatting strings
            handler = logging.StreamHandler()  # output to the console
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            MyDepth.logger_initialized = True
        else:
            self.logger = logging.getLogger(__name__)
        # arg
        self.samplename = samplename
        self.lbam_path = lbam_path
        self.abam_path = abam_path
        self.depth_path = depth_path
        self.jobs = jobs
        self.threads = threads
        self.restart = restart
        self.sequenceDict = sequenceDict


    # Longread depth
    def longread_depth(self):
        longreaddepth_log_out = ''

        # path
        longreaddepth_out =  os.path.join(self.depth_path, f'{self.samplename}.longread.depth.txt')

        # cmd
        longreaddepth_cmd = f'samtools depth -a {self.lbam_path} -@ {self.threads} -o {longreaddepth_out}'

        # Check if longread depth already exists
        if self.restart and os.path.exists(longreaddepth_out):
            self.logger.error(
                f'TGS depth file already exists: {longreaddepth_out}., skipping...')
        else:
            try:
                self.logger.error(f'Calculate TGS based alignment depth: {longreaddepth_cmd}')
                stdout, stderr, log_out = run_cmd.run(longreaddepth_cmd, env_path)
                if log_out and not os.path.exists(f'{longreaddepth_out}'):
                    longreaddepth_log_out = log_out.strip()
                self.logger.error('TGS based alignment depth calculation is done ...')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        ## depth_dict

        # Initialize the longread comparison results to build a dictionary with storage depth
        self.logger.error(f'Build TGS alignment depth index ...')
        longread_depth_dict = {key: [0] * len(value) for key, value in self.sequenceDict.items()}
        # I/O vcf
        self.longreaddepthClass = MyVCF(longreaddepth_out)
        longreaddepthOpen = self.longreaddepthClass.open_vcf_file()
        with longreaddepthOpen as f_longread:
            for line in f_longread:
                line_list = line.strip().split("\t")
                chromosome = line_list[0]
                pos = int(line_list[1])
                dep = int(line_list[2])
                longread_depth_dict[chromosome][pos - 1] = dep
        self.logger.error(f'TGS alignment depth index is done ...')

        return longread_depth_dict, longreaddepth_log_out


    # Genome depth
    def genome_depth(self):
        genomedepth_log_out = ''

        # path
        genomedepth_out = os.path.join(self.depth_path, f'{self.samplename}.genome.depth.txt')

        # cmd
        genomedepth_cmd = f'samtools depth -a {self.abam_path} -@ {self.threads} -o {genomedepth_out}'

        # Check if genome depth already exists
        if self.restart and os.path.exists(genomedepth_out):
            self.logger.error(
                f'Genome depth file already exists: {genomedepth_out}., skipping...')
        else:
            try:
                self.logger.error(f'Calculate Genome based alignment depth: {genomedepth_cmd}')
                stdout, stderr, log_out = run_cmd.run(genomedepth_cmd, env_path)
                if log_out and not os.path.exists(f'{genomedepth_out}'):
                    genomedepth_log_out = log_out.strip()
                self.logger.error('Genome based alignment depth calculation is done ...')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        ## depth_dict

        # Initialize the longread comparison results to build a dictionary with storage depth
        self.logger.error(f'Build Genome alignment depth index ...')
        genome_depth_dict = {key: [0] * len(value) for key, value in self.sequenceDict.items()}
        # I/O vcf
        self.genomedepthClass = MyVCF(genomedepth_out)
        genomedepthOpen = self.genomedepthClass.open_vcf_file()
        with genomedepthOpen as f_genome:
            for line in f_genome:
                line_list = line.strip().split("\t")
                chromosome = line_list[0]
                pos = int(line_list[1])
                dep = int(line_list[2])
                genome_depth_dict[chromosome][pos - 1] = dep
        self.logger.error(f'Genome alignment depth index is done ...')

        return genome_depth_dict, genomedepth_log_out


    # Main method with multiprocessing
    def run_depth(self):
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.jobs) as executor:
            futures = []
            depth_path_dict = {}
            if self.lbam_path:
                future_longreaddepth = executor.submit(self.longread_depth)
                futures.append(('longreaddepth', future_longreaddepth))
            if self.abam_path:
                future_genomedepth = executor.submit(self.genome_depth)
                futures.append(('genomedepth', future_genomedepth))
        for task_name, future in futures:
            depth_path, log_out = future.result()
            if log_out:
                raise SystemExit(log_out)
            depth_path_dict[task_name] = depth_path

        return depth_path_dict