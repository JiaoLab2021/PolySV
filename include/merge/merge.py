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
import shutil


# environment variables
env_path = {'PATH': os.environ.get('PATH')}
# include
include_dir = os.path.dirname(os.path.realpath(__file__))
for root, dirs, files in os.walk(os.path.join(include_dir, "include")):
    sys.path.append(root)


from make_dir import MyMkdir
from vcf_split import MyVcfsplit
from vcf_io import MyVCF
from vcf_sort import MyVcfsort


class MyMerge():
    # Flag whether global log handler has been added
    logger_initialized = False

    def __init__(self, vcf2convert_path_dict, sequenceDict, reference_path, fq_path, fa_path, samplename, jobs, threads,
                 workdir_path, restart):
        # If the log handler has not been added yet, add
        if not MyMerge.logger_initialized:
            self.logger = logging.getLogger(__name__)
            formatter = logging.Formatter('[%(asctime)s] %(message)s')  # add date and time formatting strings
            handler = logging.StreamHandler()  # output to the console
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            MyMerge.logger_initialized = True
        else:
            self.logger = logging.getLogger(__name__)

        # Value
        self.vcf2convert_path_dict = vcf2convert_path_dict
        self.sequenceDict = sequenceDict
        self.mergelarge_path = ''
        self.mergelarge_tmppath = ''

        # Parameter
        self.reference_path = reference_path
        self.fq_path = fq_path
        self.fa_path = fa_path
        self.samplename = samplename
        self.jobs = jobs
        self.threads = threads
        self.workdir_path = workdir_path
        self.restart = restart

        # Class
        self.mkdirClass = MyMkdir()


    ## Split VCF

    # Split Assemblytics2convert
    def split_Assemblytics2convert_vcf(self, mergelarge_tmppath):
        # path
        Assemblytics2convert_vcf = self.vcf2convert_path_dict['Assemblytics2convert_vcf']  # input
        Assemblytics2convertstandard_vcf = os.path.join(mergelarge_tmppath,
                                                        f'{self.samplename}.Assemblytics2convert.s.vcf')
        Assemblytics2convertlarge_vcf = os.path.join(mergelarge_tmppath,
                                                        f'{self.samplename}.Assemblytics2convert.l.vcf')

        # Check if svimasm2convert split already exists
        if self.restart and os.path.exists(Assemblytics2convertstandard_vcf) and \
                os.path.exists(Assemblytics2convertlarge_vcf):
            self.logger.error(f'The split format of Genome Assemblytics2convert vcf already exists: '
                              f'{Assemblytics2convertstandard_vcf} and {Assemblytics2convertlarge_vcf}., skipping...')
        else:
            try:
                self.logger.error(f'Start spliting the Assemblytics2convert vcf')
                # split
                self.Assemblytics2convertsplitClass = MyVcfsplit(self.threads, Assemblytics2convert_vcf,
                                                                 Assemblytics2convertstandard_vcf,
                                                                 Assemblytics2convertlarge_vcf, self.sequenceDict)
                self.Assemblytics2convertsplitClass.output()
                self.logger.error('The split format of Genome Assemblytics2convert vcf is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return Assemblytics2convertstandard_vcf, Assemblytics2convertlarge_vcf

    # Split svim-asm2convert
    def split_svimasm2convert_vcf(self, mergelarge_tmppath):
        # path
        svimasm2convert_vcf = self.vcf2convert_path_dict['svim-asm2convert_vcf']  # input
        svimasm2convertstandard_vcf = os.path.join(mergelarge_tmppath,
                                                        f'{self.samplename}.svim-asm2convert.s.vcf')
        svimasm2convertlarge_vcf = os.path.join(mergelarge_tmppath,
                                                     f'{self.samplename}.svim-asm2convert.l.vcf')

        # Check if svimasm2convert split already exists
        if self.restart and os.path.exists(svimasm2convertstandard_vcf) and \
                os.path.exists(svimasm2convertlarge_vcf):
            self.logger.error(f'The split format of Genome Assemblytics2convert vcf already exists: '
                              f'{svimasm2convertstandard_vcf} and {svimasm2convertlarge_vcf}., skipping...')
        else:
            try:
                self.logger.error(f'Start spliting the svim-asm2convert vcf')
                # split
                self.svimasm2convertsplitClass = MyVcfsplit(self.threads, svimasm2convert_vcf,
                                                                 svimasm2convertstandard_vcf,
                                                                 svimasm2convertlarge_vcf, self.sequenceDict)
                self.svimasm2convertsplitClass.output()
                self.logger.error('The split format of Genome svim-asm2convert vcf is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return svimasm2convertstandard_vcf, svimasm2convertlarge_vcf


    # Split
    def split_vcf(self):
        # Vcf dict
        standardvcf_dict = {}
        largevcf_dict = {}


        # longread
        if self.fq_path:
            longread_callers = ['sniffles', 'cuteSV', 'debreak', 'svim']
            standardvcf_dict = {longread_caller: self.vcf2convert_path_dict[f'{longread_caller}2convert_vcf']
                                for longread_caller in longread_callers}

        # path
        self.mergelarge_path = os.path.join(self.workdir_path, "04_merge")
        self.mkdirClass.make(self.mergelarge_path, self.restart)
        self.mergelarge_tmppath = os.path.join(self.mergelarge_path, "tmp_work")
        self.mkdirClass.make(self.mergelarge_tmppath, self.restart)

        # fa
        if self.fa_path:
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.jobs) as executor:
                futures_split = []

                # Assemblytics
                future_Assemblytics2convertsplit = executor.submit(self.split_Assemblytics2convert_vcf, self.mergelarge_tmppath)
                futures_split.append(('Assemblytics', future_Assemblytics2convertsplit))
                # svim-asm
                future_svimasm2convertsplit = executor.submit(self.split_svimasm2convert_vcf, self.mergelarge_tmppath)
                futures_split.append(('svim-asm', future_svimasm2convertsplit))
            # return split file path
            for task_name, future_split in futures_split:
                standard_vcf, convertlarge_vcf = future_split.result()
                standardvcf_dict[task_name] = standard_vcf
                largevcf_dict[task_name] = convertlarge_vcf

        return standardvcf_dict, largevcf_dict


    # Build configuration file
    def make_config(self):
        # Split
        standardvcf_dict, largevcf_dict = self.split_vcf()

        # Defining variables
        genomeandlongread_callers = ['svim-asm', 'sniffles', 'cuteSV', 'svim', 'debreak', 'Assemblytics']
        genome_callers = ['svim-asm', 'Assemblytics']
        longread_callers = ['sniffles', 'cuteSV', 'svim', 'debreak']
        jasminesconfig_txt = os.path.abspath(os.path.join(self.mergelarge_tmppath,
                                                          f'{self.samplename}.jasmine.s.config.txt'))
        jasminelconfig_txt = ''

        # jasminesconfig_txt
        if standardvcf_dict and largevcf_dict:
            if len(standardvcf_dict) == len(largevcf_dict):  # genome
                standardvcf_order = genome_callers
            else:
                standardvcf_order = genomeandlongread_callers
        elif standardvcf_dict:
            standardvcf_order = longread_callers
        elif largevcf_dict:
            standardvcf_order = genome_callers
        else:
            raise ValueError("No valid VCF file")

        self.logger.error("Building the Jasmine standard configuration file")
        self.jasminesconfigClass = MyVCF(jasminesconfig_txt)
        jasminesconfigWrite = self.jasminesconfigClass.write_vcf_file()
        self.logger.error(f"The jasmine standard configuration file is built: {jasminesconfig_txt}")

        with jasminesconfigWrite as f_config:
            for caller in standardvcf_order:
                if caller in standardvcf_dict:
                    f_config.write(standardvcf_dict[caller] + '\n')

        # jasminelconfig_txt
        if largevcf_dict:
            jasminelconfig_txt = os.path.abspath(os.path.join(self.mergelarge_tmppath,
                                                              f'{self.samplename}.jasmine.l.config.txt'))

            largevcf_order = genome_callers
            self.logger.error(f"Building the Jasmine large configuration file")
            self.jasminelconfigClass = MyVCF(jasminelconfig_txt)
            jasminelconfigWrite = self.jasminelconfigClass.write_vcf_file()
            self.logger.error(f"The jasmine large configuration file is built: {jasminelconfig_txt}")
            with jasminelconfigWrite as f_config:
                for caller in largevcf_order:
                    f_config.write(largevcf_dict[caller] + '\n')

        return jasminesconfig_txt, jasminelconfig_txt


    # Merge jasmine
    def jasmine_merge(self, jasmineconfig_txt, jasmineoutdir_path, jasmineoutvcf_path, support):
        jasmine_log_out = ''

        # path
        self.final_mergevcfgz_path = os.path.join(self.mergelarge_path, f'{self.samplename}.jasmine.vcf.gz')

        # cmd
        jasmine_cmd = f'jasmine file_list={jasmineconfig_txt} genome_file={self.reference_path} out_dir' \
                      f'={jasmineoutdir_path} out_file={jasmineoutvcf_path} max_dist=1000 min_support={support} ' \
                      f'threads={self.threads} --dup_to_ins --ignore_type --ignore_strand --normalize_type'
        bgzip_cmd = f'bgzip -f -@ {self.threads} {jasmineoutvcf_path} && tabix -f {jasmineoutvcf_path}.gz'

        # Check if jasmine out already exists
        if self.restart and os.path.exists(self.final_mergevcfgz_path):
            self.logger.error(
                f'Merge on jasmine already exists: {jasmineoutvcf_path}.gz '
                f'and {jasmineoutvcf_path}.gz.tbi., skipping...')
        else:
            try:
                self.logger.error(f'Start merge on jasmine: {jasmine_cmd}')
                stdout, stderr, log_out = run_cmd.run(jasmine_cmd, env_path)
                if log_out and not os.path.exists(f'{jasmineoutvcf_path}'):
                    jasmine_log_out = log_out.strip()
                self.logger.error('Merge on jasmine is done.')

                # sort
                self.jasminesortClass = MyVcfsort(jasmineoutvcf_path)
                self.jasminesortClass.sort()

                # format(For bcftools)
                jasmineClass = MyVCF(jasmineoutvcf_path)  # input
                jasmineOpen = jasmineClass.open_vcf_file()
                jasmineformat_path = os.path.abspath(f'{jasmineoutvcf_path}.format')  # output
                jasmineformatClass = MyVCF(jasmineformat_path)
                jasmineformatWrite = jasmineformatClass.write_vcf_file()
                with jasmineOpen as f_in:
                    for line in f_in:
                        if "#" in line:
                            jasmineformatWrite.write(line)
                        else:
                            info = line.strip()
                            infoList = info.split("\t")
                            INFOCol = infoList[7]

                            INFOCol_feature = {"SUPP=", "END=", "SVLEN=", "OLDTYPE=", "IDLIST="}
                            INFOCol_list = [
                                item for item in INFOCol.split(';') if
                                any(keyword in item for keyword in INFOCol_feature)
                            ]
                            infoList[7] = ';'.join(INFOCol_list)
                            jasmineformatWrite.write("\t".join(infoList) + "\n")
                if jasmineformat_path:
                    os.rename(jasmineformat_path, jasmineoutvcf_path)

                # bgzip
                self.logger.error(f'Compress and build index: {bgzip_cmd}')
                run_cmd.run(bgzip_cmd, env_path)

            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return jasmine_log_out


    # Main method with multiprocessing
    def run_merge(self):
        # Make_config
        jasminesconfig_txt, jasminelconfig_txt = self.make_config()
        standard_support = 2
        large_support = 1

        # Sort
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.jobs) as executor:

            # standard
            jasminesoutdir_path = os.path.join(self.mergelarge_tmppath, f'{self.samplename}.jasmine.s.out')
            jasminesoutvcf_path = os.path.join(self.mergelarge_path, f'{self.samplename}.jasmine.s.vcf')
            future_jasmines = executor.submit(self.jasmine_merge, jasminesconfig_txt,
                                              jasminesoutdir_path, jasminesoutvcf_path, standard_support)
            log_jasminesout = future_jasmines.result()
            if log_jasminesout:
                raise SystemExit(log_jasminesout)


            # Whether large
            if jasminesconfig_txt and jasminelconfig_txt:
                jasmineloutdir_path = os.path.join(self.mergelarge_tmppath, f'{self.samplename}.jasmine.l.out')
                jasmineloutvcf_path = os.path.join(self.mergelarge_path, f'{self.samplename}.jasmine.l.vcf')
                future_jasminel = executor.submit(self.jasmine_merge, jasminelconfig_txt,
                                                  jasmineloutdir_path, jasmineloutvcf_path, large_support)
                log_jasminelout = future_jasminel.result()
                if log_jasminelout:
                    raise SystemExit(log_jasminelout)

                # bcftools
                bcftools_cmd = f'bcftools concat --threads {self.threads} -a {jasminesoutvcf_path}.gz ' \
                               f'{jasmineloutvcf_path}.gz | bcftools norm -m +any - -Oz -o {self.final_mergevcfgz_path}'
                if self.restart and os.path.exists(self.final_mergevcfgz_path):
                    self.logger.error(f'Final merged VCF file already exists: {self.final_mergevcfgz_path}, skipping...')
                else:
                    try:
                        self.logger.error(f"Generate the final merged VCF file: {bcftools_cmd}")
                        stdout, stderr, log_out = run_cmd.run(bcftools_cmd, env_path)
                        if log_out and not os.path.exists(self.final_mergevcfgz_path):
                            raise SystemExit(log_out)
                        self.logger.error('Final merged VCF file is done.')
                    except Exception as e:
                        self.logger.error(f"An error occurred: {e}")

                remove_files = [f'{jasminesoutvcf_path}.gz', f'{jasminesoutvcf_path}.gz.tbi',
                                f'{jasmineloutvcf_path}.gz', f'{jasmineloutvcf_path}.gz.tbi']
            else:
                if self.restart and os.path.exists(self.final_mergevcfgz_path):
                    self.logger.error(f'Final merged VCF file already exists: {self.final_mergevcfgz_path}, skipping...')
                else:
                    os.rename(f'{jasminesoutvcf_path}.gz', self.final_mergevcfgz_path)

                remove_files = [f'{jasminesoutvcf_path}.gz', f'{jasminesoutvcf_path}.gz.tbi']

            # remove
            if os.path.exists(self.final_mergevcfgz_path):
                for remove_file in remove_files:
                    if os.path.exists(remove_file):
                        os.remove(remove_file)
                shutil.rmtree(self.mergelarge_tmppath)

        return self.final_mergevcfgz_path