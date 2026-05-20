#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"


import os
import sys
import logging
import concurrent.futures
import traceback


# environment variables
env_path = {'PATH': os.environ.get('PATH')}
# include
include_dir = os.path.dirname(os.path.realpath(__file__))
for root, dirs, files in os.walk(os.path.join(include_dir, "include")):
    sys.path.append(root)

from make_dir import MyMkdir
from sniffles2vcf import MySniffles2convert
from cuteSV2vcf import MyCuteSV2convert
from debreak2vcf import MyDebreak2convert
from svim2vcf import MySvim2convert
from Assemblytics2vcf import MyAssemblytics2convert
from svimasm2vcf import MySvimasm2convert


# get parser
class MyVcf2Convert:
    # Flag whether global log handler has been added
    logger_initialized = False

    def __init__(self, call_path_dict, sequenceDict, reference_path, fq_path, fa_path, fqdt, samplename, jobs, threads,
                 workdir_path, restart):
        # If the log handler has not been added yet, add
        if not MyVcf2Convert.logger_initialized:
            self.logger = logging.getLogger(__name__)
            formatter = logging.Formatter('[%(asctime)s] %(message)s')  # add date and time formatting strings
            handler = logging.StreamHandler()  # output to the console
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            MyVcf2Convert.logger_initialized = True
        else:
            self.logger = logging.getLogger(__name__)

        # Value
        self.call_path_dict = call_path_dict
        self.sequenceDict = sequenceDict

        # Parameter
        self.reference_path = reference_path
        self.fq_path = fq_path
        self.fa_path = fa_path
        self.fqdt = fqdt
        self.samplename = samplename
        self.jobs = jobs
        self.threads = threads
        self.workdir_path = workdir_path
        self.restart = restart

        # Class
        self.mkdirClass = MyMkdir()


    # Convert sniffles
    def convert_sniffles_longread(self, longreadsniffles_prefix, vcf2convertlongread_path):
        # path
        vcf2convertlongreadsniffles_path = os.path.join(vcf2convertlongread_path, "sniffles")
        self.mkdirClass.make(vcf2convertlongreadsniffles_path, self.restart)
        sniffles2convert_prefix = os.path.join(vcf2convertlongreadsniffles_path, self.samplename)

        sniffles_vcfgz = f'{longreadsniffles_prefix}.sniffles.vcf.gz'
        sniffles2convert_vcf = os.path.abspath(f'{sniffles2convert_prefix}.sniffles2convert.vcf')

        # Check if sniffles2convert out already exists
        if self.restart and os.path.exists(sniffles2convert_vcf):
            self.logger.error(f'The convert format of TGS sniffles already exists: '
                              f'{sniffles2convert_vcf}., skipping...')
        else:
            try:
                self.logger.error(f'Start converting the output file format of sniffles')
                # sniffles2vcf
                self.sniffles2convertClass = MySniffles2convert(self.reference_path, self.samplename, self.threads,
                                                                sniffles_vcfgz, sniffles2convert_vcf,
                                                                self.sequenceDict)
                self.sniffles2convertClass.output()
                self.logger.error('The convert format of TGS sniffles is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return sniffles2convert_vcf


    # Convert cuteSV
    def convert_cuteSV_longread(self, longreadcuteSV_prefix, vcf2convertlongread_path):
        # path
        vcf2convertlongreadcuteSV_path = os.path.join(vcf2convertlongread_path, "cuteSV")
        self.mkdirClass.make(vcf2convertlongreadcuteSV_path, self.restart)
        cuteSV2convert_prefix = os.path.join(vcf2convertlongreadcuteSV_path, self.samplename)

        cuteSV_vcfgz = f'{longreadcuteSV_prefix}.cuteSV.vcf.gz'
        cuteSV2convert_vcf = os.path.abspath(f'{cuteSV2convert_prefix}.cuteSV2convert.vcf')

        # Check if cuteSV2convert out already exists
        if self.restart and os.path.exists(cuteSV2convert_vcf):
            self.logger.error(f'The convert format of TGS cuteSV already exists: '
                              f'{cuteSV2convert_vcf}., skipping...')
        else:
            try:
                self.logger.error(f'Start converting the output file format of cuteSV')
                # cuteSV2vcf
                self.cuteSV2convertClass = MyCuteSV2convert(self.reference_path, self.samplename, self.threads,
                                                                cuteSV_vcfgz, cuteSV2convert_vcf,
                                                                self.sequenceDict)
                self.cuteSV2convertClass.output()
                self.logger.error('The convert format of TGS cuteSV is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")
                self.logger.error(traceback.format_exc())

        return cuteSV2convert_vcf


    # Convert debreak
    def convert_debreak_longread(self, longreaddebreak_prefix, vcf2convertlongread_path):
        # path
        vcf2convertlongreaddebreak_path = os.path.join(vcf2convertlongread_path, "debreak")
        self.mkdirClass.make(vcf2convertlongreaddebreak_path, self.restart)
        debreak2convert_prefix = os.path.join(vcf2convertlongreaddebreak_path, self.samplename)

        debreak_vcfgz = f'{longreaddebreak_prefix}.debreak.vcf.gz'
        debreak2convert_vcf = os.path.abspath(f'{debreak2convert_prefix}.debreak2convert.vcf')

        # Check if debreak2convert out already exists
        if self.restart and os.path.exists(debreak2convert_vcf):
            self.logger.error(f'The convert format of TGS debreak already exists: '
                              f'{debreak2convert_vcf}., skipping...')
        else:
            try:
                self.logger.error(f'Start converting the output file format of debreak')
                # debreak2vcf
                self.debreak2convertClass = MyDebreak2convert(self.reference_path, self.samplename, self.threads,
                                                            debreak_vcfgz, debreak2convert_vcf,
                                                            self.sequenceDict)
                self.debreak2convertClass.output()
                self.logger.error('The convert format of TGS debreak is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return debreak2convert_vcf


    # Convert svim
    def convert_svim_longread(self, longreadsvim_prefix, vcf2convertlongread_path):
        # path
        vcf2convertlongreadsvim_path = os.path.join(vcf2convertlongread_path, "svim")
        self.mkdirClass.make(vcf2convertlongreadsvim_path, self.restart)
        svim2convert_prefix = os.path.join(vcf2convertlongreadsvim_path, self.samplename)

        svim_vcfgz = f'{longreadsvim_prefix}.svim.vcf.gz'
        svim2convert_vcf = os.path.abspath(f'{svim2convert_prefix}.svim2convert.vcf')

        # Check if svim2convert out already exists
        if self.restart and os.path.exists(svim2convert_vcf):
            self.logger.error(f'The convert format of TGS svim already exists: '
                              f'{svim2convert_vcf}., skipping...')
        else:
            try:
                self.logger.error(f'Start converting the output file format of svim')
                # svim2vcf
                self.svim2convertClass = MySvim2convert(self.reference_path, self.samplename, self.threads,
                                                              svim_vcfgz, svim2convert_vcf,
                                                              self.sequenceDict)
                self.svim2convertClass.output()
                self.logger.error('The convert format of TGS svim is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return svim2convert_vcf


    # Convert Assemblytics
    def convert_Assemblytics_genome(self, genomeAssemblytics_prefix, vcf2convertgenome_path):
        # path
        vcf2convertgenomeAssemblytics_path = os.path.join(vcf2convertgenome_path, "Assemblytics")
        self.mkdirClass.make(vcf2convertgenomeAssemblytics_path, self.restart)
        Assemblytics2convert_prefix = os.path.join(vcf2convertgenomeAssemblytics_path, self.samplename)

        Assemblytics_bedgz = f'{genomeAssemblytics_prefix}.Assemblytics.bed.gz'
        Assemblytics2convert_vcf = os.path.abspath(f'{Assemblytics2convert_prefix}.Assemblytics2convert.vcf')

        # Check if Assemblytics2convert out already exists
        if self.restart and os.path.exists(Assemblytics2convert_vcf):
            self.logger.error(f'The convert format of TGS Assemblytics already exists: '
                              f'{Assemblytics2convert_vcf}., skipping...')
        else:
            try:
                self.logger.error(f'Start converting the output file format of Assemblytics')
                # Assemblytics2vcf
                self.Assemblytics2convertClass = MyAssemblytics2convert(self.reference_path, self.fa_path,
                                                                        self.samplename, self.threads,
                                                                        Assemblytics_bedgz, Assemblytics2convert_vcf,
                                                                        self.sequenceDict)
                self.Assemblytics2convertClass.output()
                self.logger.error('The convert format of Genome Assemblytics is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return Assemblytics2convert_vcf


    # Convert svim-asm
    def convert_svimasm_genome(self, genomesvimasm_prefix, vcf2convertgenome_path):
        # path
        vcf2convertgenomesvimasm_path = os.path.join(vcf2convertgenome_path, "svim-asm")
        self.mkdirClass.make(vcf2convertgenomesvimasm_path, self.restart)
        svimasm2convert_prefix = os.path.join(vcf2convertgenomesvimasm_path, self.samplename)

        svimasm_vcfgz = f'{genomesvimasm_prefix}.svim-asm.vcf.gz'
        svimasm2convert_vcf = os.path.abspath(f'{svimasm2convert_prefix}.svim-asm2convert.vcf')

        # Check if svimasm2convert out already exists
        if self.restart and os.path.exists(svimasm2convert_vcf):
            self.logger.error(f'The convert format of Genome svim-asm already exists: '
                              f'{svimasm2convert_vcf}., skipping...')
        else:
            try:
                self.logger.error(f'Start converting the output file format of svim-asm')
                # svimasm2vcf
                self.svimasm2convertClass = MySvimasm2convert(self.reference_path, self.samplename, self.threads,
                                                        svimasm_vcfgz, svimasm2convert_vcf,
                                                        self.sequenceDict)
                self.svimasm2convertClass.output()
                self.logger.error('The convert format of TGS svimasm is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return svimasm2convert_vcf


    # Main method with multiprocessing
    def run_vcf2convert(self):
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.jobs) as executor:
            futures = []
            vcf2convert_path_dict = {}
            if self.fq_path:
                vcf2convertlongread_path = os.path.join(self.workdir_path, "03_vcf2convert", "01_longread")
                self.mkdirClass.make(vcf2convertlongread_path, self.restart)
                # sniffles
                longreadsniffles_prefix = self.call_path_dict['longreadsniffles_prefix']
                future_longreadsniffles2convert = executor.submit(self.convert_sniffles_longread,
                                                                  longreadsniffles_prefix, vcf2convertlongread_path)
                futures.append(('sniffles2convert_vcf', future_longreadsniffles2convert))
                # cuteSV
                longreadcuteSV_prefix = self.call_path_dict['longreadcuteSV_prefix']
                future_longreadcuteSV2convert = executor.submit(self.convert_cuteSV_longread,
                                                                  longreadcuteSV_prefix, vcf2convertlongread_path)
                futures.append(('cuteSV2convert_vcf', future_longreadcuteSV2convert))
                # debreak
                longreaddebreak_prefix = self.call_path_dict['longreaddebreak_prefix']
                future_longreaddebreak2convert = executor.submit(self.convert_debreak_longread,
                                                                longreaddebreak_prefix, vcf2convertlongread_path)
                futures.append(('debreak2convert_vcf', future_longreaddebreak2convert))
                # svim
                longreadsvim_prefix = self.call_path_dict['longreadsvim_prefix']
                future_longreadsvim2convert = executor.submit(self.convert_svim_longread,
                                                                 longreadsvim_prefix, vcf2convertlongread_path)
                futures.append(('svim2convert_vcf', future_longreadsvim2convert))
            if self.fa_path:
                vcf2convertgenome_path = os.path.join(self.workdir_path, "03_vcf2convert", "02_genome")
                self.mkdirClass.make(vcf2convertgenome_path, self.restart)
                # Assemblytics
                genomeAssemblytics_prefix = self.call_path_dict['genomeAssemblytics_prefix']
                future_genomeAssemblytics2convert = executor.submit(self.convert_Assemblytics_genome,
                                                                  genomeAssemblytics_prefix, vcf2convertgenome_path)
                futures.append(('Assemblytics2convert_vcf', future_genomeAssemblytics2convert))
                # svim-asm
                genomesvimasm_prefix = self.call_path_dict['genomesvimasm_prefix']
                future_genomesvimasm2convert = executor.submit(self.convert_svimasm_genome,
                                                                    genomesvimasm_prefix, vcf2convertgenome_path)
                futures.append(('svim-asm2convert_vcf', future_genomesvimasm2convert))

        for task_name, future in futures:
            vcf2convert_path = future.result()
            vcf2convert_path_dict[task_name] = vcf2convert_path

        return vcf2convert_path_dict