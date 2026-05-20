#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"


import os
import sys
import logging
import run_cmd
import concurrent.futures


# environment variables
env_path = {'PATH': os.environ.get('PATH')}
# include
include_dir = os.path.dirname(os.path.realpath(__file__))
for root, dirs, files in os.walk(os.path.join(include_dir, "include")):
    sys.path.append(root)


from make_dir import MyMkdir


class MyAlign():
    # Flag whether global log handler has been added
    logger_initialized = False

    def __init__(self, reference_path, fq_path, fa_path, fqdt, samplename, jobs, threads, workdir_path, restart):
        # If the log handler has not been added yet, add
        if not MyAlign.logger_initialized:
            self.logger = logging.getLogger(__name__)
            formatter = logging.Formatter('[%(asctime)s] %(message)s')  # add date and time formatting strings
            handler = logging.StreamHandler()  # output to the console
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            MyAlign.logger_initialized = True
        else:
            self.logger = logging.getLogger(__name__)

        # Parameter
        self.reference_path = reference_path
        self.fq_path = fq_path
        self.fa_path = fa_path
        self.fqdt = fqdt
        self.jobs = jobs
        self.threads = threads
        self.samplename = samplename
        self.workdir_path = workdir_path
        self.restart = restart

        # Class
        self.mkdirClass = MyMkdir()


    # Long-read alignment
    def run_minimap2_longread(self, alignlongread_path):
        longreadminimap2_log_out = ''

        minimap2_fqdt_dic = {'ccs': 'map-hifi', 'clr': 'map-pb', 'ont': 'map-ont'}
        minimap2_fqdt = minimap2_fqdt_dic.get(self.fqdt)
        if not minimap2_fqdt:
            raise ValueError(
                f"Invalid fqdt value: {self.fqdt}. Expected one of: {', '.join(minimap2_fqdt_dic.keys())}")

        # path
        alignlongreadminimap2_path = os.path.join(alignlongread_path, "minimap2")
        self.mkdirClass.make(alignlongreadminimap2_path, self.restart)
        longreadminimap2_prefix = os.path.join(alignlongreadminimap2_path, self.samplename)
        alignlongreadminimap2_bam = os.path.join(f'{longreadminimap2_prefix}.bam')

        # cmd
        longreadminimap2_cmd = f'minimap2 -t {self.threads} -ax {minimap2_fqdt} {self.reference_path} ' \
                               f'{self.fq_path} | samtools sort -@ {self.threads} --write-index ' \
                               f'-o {alignlongreadminimap2_bam}'

        # Check if alignment already exists
        if self.restart and os.path.exists(alignlongreadminimap2_bam):
            self.logger.error(f'TGS alignment on minimap2 already exists: {alignlongreadminimap2_bam}, skipping...')
        else:
            try:
                self.logger.error(f'Start TGS alignment on minimap2: {longreadminimap2_cmd}')
                stdout, stderr, log_out = run_cmd.run(longreadminimap2_cmd, env_path)
                if log_out:
                    longreadminimap2_log_out = log_out.strip()
                self.logger.error('TGS alignment on minimap2 is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return longreadminimap2_prefix, longreadminimap2_log_out


    # Genome minimap2 alignment
    def run_minimap2_genome(self, aligngenome_path):
        genomeminimap2_log_out = ''

        # path
        aligngenomeminimap2_path = os.path.join(aligngenome_path, "minimap2")
        self.mkdirClass.make(aligngenomeminimap2_path, self.restart)
        genomeminimap2_prefix = os.path.join(aligngenomeminimap2_path, self.samplename)
        aligngenomeminimap2_bam = os.path.join(f'{genomeminimap2_prefix}.bam')

        # cmd
        genomeminimap2_cmd = f'minimap2 -t {self.threads} -ax asm5 --cs -r2k {self.reference_path} ' \
                             f'{self.fa_path} | samtools sort -@ {self.threads} --write-index ' \
                             f'-o {aligngenomeminimap2_bam}'

        # Check if alignment already exists
        if self.restart and os.path.exists(aligngenomeminimap2_bam):
            self.logger.error(
                f'Genome alignment on minimap2 already exists: {aligngenomeminimap2_bam}, skipping...')
        else:
            try:
                self.logger.error(f'Start Genome alignment on minimap2: {genomeminimap2_cmd}')
                stdout, stderr, log_out = run_cmd.run(genomeminimap2_cmd, env_path)
                if log_out:
                    genomeminimap2_log_out=log_out.strip()
                self.logger.error('Genome alignment on minimap2 is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return genomeminimap2_prefix, genomeminimap2_log_out


    # Genome nucmer alignment
    def run_nucmer_genome(self, aligngenome_path):
        genomenucmer_log_out = ''

        # path
        aligngenomenucmer_path = os.path.join(aligngenome_path, "nucmer")
        self.mkdirClass.make(aligngenomenucmer_path, self.restart)

        genomenucmer_prefix = os.path.join(aligngenomenucmer_path, self.samplename)
        aligngenomenucmer_delta = os.path.abspath(f'{genomenucmer_prefix}.delta')
        aligngenomenucmer_filterdelta = os.path.abspath(f'{genomenucmer_prefix}.filter.delta')
        aligngenomenucmer_coords = os.path.abspath(f'{genomenucmer_prefix}.coords')

        # cmd
        genomenucmer_cmd = f'nucmer -t {self.threads} --batch 1 -c 500 -b 500 -l 100 {self.reference_path} ' \
                           f'{self.fa_path} --prefix={genomenucmer_prefix} && delta-filter -i 90 -l 100 ' \
                           f'{aligngenomenucmer_delta} 1>{aligngenomenucmer_filterdelta} && show-coords ' \
                           f'-THrd {aligngenomenucmer_filterdelta} 1>{aligngenomenucmer_coords}'

        # Check if alignment already exists
        if self.restart and os.path.exists(aligngenomenucmer_delta) and os.path.exists\
                    (aligngenomenucmer_filterdelta) and os.path.exists(aligngenomenucmer_coords):
            self.logger.error(f'Genome alignment on nucmer already exists: {aligngenomenucmer_coords}, skipping...')
        else:
            try:
                self.logger.error(f'Start Genome alignment on nucmer: {genomenucmer_cmd}')
                stdout, stderr, log_out = run_cmd.run(genomenucmer_cmd, env_path)
                if log_out:
                    genomenucmer_log_out = log_out.strip()
                self.logger.error('Genome alignment on nucmer is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return genomenucmer_prefix, genomenucmer_log_out


    # Main method with multiprocessing
    def run_align(self):
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.jobs) as executor:
            futures = []
            align_path_dict = {}
            if self.fq_path:
                alignlongread_path = os.path.abspath(os.path.join(self.workdir_path, "01_align", "01_longread"))
                self.mkdirClass.make(alignlongread_path, self.restart)

                future_longreadminimap2 = executor.submit(self.run_minimap2_longread, alignlongread_path)
                futures.append(('longreadminimap2_prefix', future_longreadminimap2))

            if self.fa_path:
                aligngenome_path = os.path.abspath(os.path.join(self.workdir_path, "01_align", "02_genome"))
                self.mkdirClass.make(aligngenome_path, self.restart)
                future_genomeminimap2 = executor.submit(self.run_minimap2_genome, aligngenome_path)
                futures.append(('genomeminimap2_prefix', future_genomeminimap2))

                future_genomenucmer = executor.submit(self.run_nucmer_genome, aligngenome_path)
                futures.append(('genomenucmer_prefix', future_genomenucmer))

            for task_name, future in futures:
                align_path, log_out = future.result()
                if log_out:
                    raise SystemExit(log_out)
                align_path_dict[task_name] = align_path

        return align_path_dict