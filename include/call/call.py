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
import shutil

class MyCall():
    # Flag whether global log handler has been added
    logger_initialized = False

    def __init__(self, align_path_dict, reference_path, fq_path, fa_path, fqdt, samplename, jobs, threads,
                 workdir_path, restart):
        # If the log handler has not been added yet, add
        if not MyCall.logger_initialized:
            self.logger = logging.getLogger(__name__)
            formatter = logging.Formatter('[%(asctime)s] %(message)s')  # add date and time formatting strings
            handler = logging.StreamHandler()  # output to the console
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            MyCall.logger_initialized = True
        else:
            self.logger = logging.getLogger(__name__)

        # Value
        self.align_path_dict = align_path_dict

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


    # Long-read call on sniffles
    def run_sniffles_longread(self, longreadminimap2_prefix, calllongread_path):
        sniffles_log_out = ''

        # path
        calllongreadsniffles_path = os.path.join(calllongread_path, "sniffles")
        self.mkdirClass.make(calllongreadsniffles_path, self.restart)
        longreadsniffles_prefix = os.path.join(calllongreadsniffles_path, self.samplename)

        calllongreadsniffles_vcf = os.path.abspath(f'{longreadsniffles_prefix}.sniffles.vcf')
        calllongreadsniffles_vcfgz = os.path.abspath(f'{longreadsniffles_prefix}.sniffles.vcf.gz')

        # cmd
        sniffles_cmd = f'sniffles -t {self.threads} -i {longreadminimap2_prefix}.bam -v {calllongreadsniffles_vcf}'
        bgzip_cmd = f'bgzip -f -@ {self.threads} {calllongreadsniffles_vcf}'
        longreadsniffles_cmd = f'{sniffles_cmd} && {bgzip_cmd}'

        # Check if sniffles out already exists
        if self.restart and os.path.exists(calllongreadsniffles_vcfgz):
            self.logger.error(f'TGS callset on sniffles already exists: {calllongreadsniffles_vcfgz}, skipping...')
        else:
            try:
                self.logger.error(f'Start TGS callset on sniffles: {longreadsniffles_cmd}')
                stdout, stderr, log_out = run_cmd.run(longreadsniffles_cmd, env_path)
                if log_out:
                    sniffles_log_out = log_out.strip()
                else:
                    self.logger.error('TGS callset on sniffles is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return longreadsniffles_prefix, sniffles_log_out


    # Long-read call on cuteSV
    def run_cuteSV_longread(self, longreadminimap2_prefix, calllongread_path):
        cuteSV_log_out = ''

        # path
        calllongreadcuteSV_path = os.path.join(calllongread_path, "cuteSV")
        self.mkdirClass.make(calllongreadcuteSV_path, self.restart)
        longreadcuteSV_prefix = os.path.join(calllongreadcuteSV_path, self.samplename)

        calllongreadcuteSV_vcf = os.path.abspath(f'{longreadcuteSV_prefix}.cuteSV.vcf')
        calllongreadcuteSV_vcfgz = os.path.abspath(f'{longreadcuteSV_prefix}.cuteSV.vcf.gz')

        # cmd
        common_cmd = f'cuteSV -t {self.threads} ' \
                     f'--min_support 5 --min_size 50 --max_size 10000000 --genotype {longreadminimap2_prefix}.bam ' \
                     f'{self.reference_path} {calllongreadcuteSV_vcf} {calllongreadcuteSV_path}'
        bgzip_cmd = f'bgzip -f -@ {self.threads} {calllongreadcuteSV_vcf}'
        if self.fqdt == 'ccs':
            specific_cmd = '--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 ' \
                           '--max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5'
        elif self.fqdt == 'clr':
            specific_cmd = '--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 ' \
                           '--max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5'
        elif self.fqdt == 'ont':
            specific_cmd = '--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 ' \
                           '--max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3'
        else:
            raise ValueError(f"Unknown fqdt type: {self.fqdt}")
        longreadcuteSV_cmd = f'{common_cmd} {specific_cmd} && {bgzip_cmd}'

        # Check if cuteSV out already exists
        if self.restart and os.path.exists(calllongreadcuteSV_vcfgz):
            self.logger.error(f'TGS callset on cuteSV already exists: {calllongreadcuteSV_vcfgz}., skipping...')
        else:
            try:
                self.logger.error(f'Start TGS callset on cuteSV: {longreadcuteSV_cmd}')
                stdout, stderr, log_out = run_cmd.run(longreadcuteSV_cmd, env_path)
                if log_out:
                    cuteSV_log_out = log_out.strip()
                self.logger.error('TGS callset on cuteSV is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return longreadcuteSV_prefix, cuteSV_log_out


    # Long-read call on debreak
    def run_debreak_longread(self, longreadminimap2_prefix, calllongread_path):
        debreak_log_out = ''

        # path
        calllongreaddebreak_path = os.path.join(calllongread_path, "debreak")
        self.mkdirClass.make(calllongreaddebreak_path, self.restart)
        longreaddebreak_prefix = os.path.join(calllongreaddebreak_path, self.samplename)
        calllongreaddebreaktmp_path = os.path.join(calllongreaddebreak_path, "tmp_work")

        calllongreaddebreak_vcf = os.path.join(calllongreaddebreaktmp_path, f'debreak.vcf')
        calllongreaddebreak_vcfgz = os.path.abspath(f'{longreaddebreak_prefix}.debreak.vcf.gz')

        # cmd
        debreak_cmd = f'debreak --thread {self.threads} --bam {longreadminimap2_prefix}.bam --ref ' \
                      f'{self.reference_path} --min_size 50 --max_size 10000000 --min_support 5 --outpath ' \
                      f'{calllongreaddebreaktmp_path} --rescue_large_ins --rescue_dup --poa'
        bgzip_cmd = f'bgzip -f -@ {self.threads} {calllongreaddebreak_vcf}'
        longreaddebreak_cmd = f'{debreak_cmd} && {bgzip_cmd}'

        # Check if debreak out already exists
        if self.restart and os.path.exists(calllongreaddebreak_vcfgz):
            self.logger.error(f'TGS callset on debreak already exists: {calllongreaddebreak_vcfgz}., skipping...')
        else:
            try:
                self.mkdirClass.make(calllongreaddebreaktmp_path, self.restart)
                self.logger.error(f'Start TGS callset on debreak: {longreaddebreak_cmd}')
                stdout, stderr, log_out = run_cmd.run(longreaddebreak_cmd, env_path)
                if log_out:
                    debreak_log_out = log_out.strip()
                if f'{calllongreaddebreak_vcf}.gz':
                    os.rename(f'{calllongreaddebreak_vcf}.gz', f'{calllongreaddebreak_vcfgz}')
                    shutil.rmtree(calllongreaddebreaktmp_path)
                self.logger.error('TGS callset on debreak is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return longreaddebreak_prefix, debreak_log_out


    # Long-read call on svim
    def run_svim_longread(self, longreadminimap2_prefix, calllongread_path):
        svim_log_out = ''

        # path
        calllongreadsvim_path = os.path.join(calllongread_path, "svim")
        self.mkdirClass.make(calllongreadsvim_path, self.restart)
        longreadsvim_prefix = os.path.join(calllongreadsvim_path, self.samplename)
        calllongreadsvimtmp_path = os.path.join(calllongreadsvim_path, "tmp_work")

        calllongreadsvim_vcf = os.path.join(calllongreadsvimtmp_path, f'variants.vcf')
        calllongreadsvim_vcfgz = os.path.abspath(f'{longreadsvim_prefix}.svim.vcf.gz')

        # cmd
        svim_cmd = f'svim alignment --min_sv_size 50 --max_sv_size 10000000 {calllongreadsvimtmp_path} ' \
                   f'{longreadminimap2_prefix}.bam {self.reference_path}'
        bgzip_cmd = f'bgzip -f -@ {self.threads} {calllongreadsvim_vcf}'
        longreadsvim_cmd = f'{svim_cmd} && {bgzip_cmd}'

        # Check if svim out already exists
        if self.restart and os.path.exists(calllongreadsvim_vcfgz):
            self.logger.error(f'TGS callset on svim already exists: {calllongreadsvim_vcfgz}., skipping...')
        else:
            try:
                self.mkdirClass.make(calllongreadsvimtmp_path, self.restart)
                self.logger.error(f'Start TGS callset on svim: {longreadsvim_cmd}')
                stdout, stderr, log_out = run_cmd.run(longreadsvim_cmd, env_path)
                if log_out and not os.path.exists(f'{calllongreadsvim_vcf}.gz'):
                    svim_log_out = log_out.strip()
                if f'{calllongreadsvim_vcf}.gz':
                    os.rename(f'{calllongreadsvim_vcf}.gz', f'{calllongreadsvim_vcfgz}')
                    shutil.rmtree(calllongreadsvimtmp_path)
                self.logger.error('TGS callset on svim is done.')
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return longreadsvim_prefix, svim_log_out


    # Assemblytics
    def run_Assemblytics_genome(self, genomenucmer_prefix, callgenome_path):
        Assemblytics_log_out = ''

        # path
        callgenomeAssemblytics_path = os.path.join(callgenome_path, "Assemblytics")
        self.mkdirClass.make(callgenomeAssemblytics_path, self.restart)
        genomeAssemblytics_prefix = os.path.join(callgenomeAssemblytics_path, self.samplename)
        callgenomeAssemblyticstmp_path = os.path.join(callgenomeAssemblytics_path, "tmp_work")
        callgenomeAssemblyticstmpsample_path = os.path.join(callgenomeAssemblyticstmp_path, self.samplename)

        callgenomeAssemblytics_bed = os.path.join(callgenomeAssemblyticstmp_path,
                                                  f'{self.samplename}.Assemblytics_structural_variants.bed')
        callgenomeAssemblytics_bedgz = os.path.abspath(f'{genomeAssemblytics_prefix}.Assemblytics.bed.gz')

        # cmd
        Assemblytics_cmd = f'Assemblytics {genomenucmer_prefix}.filter.delta {callgenomeAssemblyticstmpsample_path} ' \
                           f'10000 50 10000000'
        bgzip_cmd = f'bgzip -f -@ {self.threads} {callgenomeAssemblytics_bed}'
        genomeAssemblytics_cmd = f'{Assemblytics_cmd} && {bgzip_cmd}'

        # Check if Assemblytics out already exists
        if self.restart and os.path.exists(callgenomeAssemblytics_bedgz):
            self.logger.error(f'TGS callset on Assemblytics already exists: {callgenomeAssemblytics_bedgz}., skipping...')
        else:
            try:
                self.mkdirClass.make(callgenomeAssemblyticstmp_path, self.restart)
                self.logger.error(f'Start TGS callset on Assemblytics: {genomeAssemblytics_cmd}')
                stdout, stderr, log_out = run_cmd.run(genomeAssemblytics_cmd, env_path)
                if log_out and not os.path.exists(f'{callgenomeAssemblytics_bed}.gz'):
                    Assemblytics_log_out = log_out.strip()
                self.logger.error('TGS callset on Assemblytics is done.')
                if f'{callgenomeAssemblytics_bed}.gz':
                    os.rename(f'{callgenomeAssemblytics_bed}.gz', f'{callgenomeAssemblytics_bedgz}')
                    shutil.rmtree(callgenomeAssemblyticstmp_path)
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return genomeAssemblytics_prefix, Assemblytics_log_out


    # svimasm
    def run_svimasm_genome(self, genomeminimap2_prefix, callgenome_path):
        svimasm_log_out = ''

        # path
        callgenomesvimasm_path = os.path.join(callgenome_path, "svim-asm")
        self.mkdirClass.make(callgenomesvimasm_path, self.restart)
        genomesvimasm_prefix = os.path.join(callgenomesvimasm_path, self.samplename)
        callgenomesvimasmtmp_path = os.path.join(callgenomesvimasm_path, "tmp_work")

        calllongreadsvimasm_vcf = os.path.join(callgenomesvimasmtmp_path, f'variants.vcf')
        calllongreadsvimasm_vcfgz = os.path.abspath(f'{genomesvimasm_prefix}.svim-asm.vcf.gz')

        # cmd
        svimasm_cmd = f'svim-asm haploid {callgenomesvimasmtmp_path} {genomeminimap2_prefix}.bam {self.reference_path} ' \
                   f'--min_mapq 30 --min_sv_size 50 --max_sv_size 10000000'
        bgzip_cmd = f'bgzip -f -@ {self.threads} {calllongreadsvimasm_vcf}'
        genomesvimasm_cmd = f'{svimasm_cmd} && {bgzip_cmd}'

        # Check if svimasm out already exists
        if self.restart and os.path.exists(calllongreadsvimasm_vcfgz):
            self.logger.error(
                f'TGS callset on svimasm already exists: {calllongreadsvimasm_vcfgz}., skipping...')
        else:
            try:
                self.mkdirClass.make(callgenomesvimasmtmp_path, self.restart)
                self.logger.error(f'Start TGS callset on svimasm: {genomesvimasm_cmd}')
                stdout, stderr, log_out = run_cmd.run(genomesvimasm_cmd, env_path)
                if log_out and not os.path.exists(f'{calllongreadsvimasm_vcf}.gz'):
                    svimasm_log_out = log_out.strip()
                self.logger.error('TGS callset on svimasm is done.')
                if f'{calllongreadsvimasm_vcf}.gz':
                    os.rename(f'{calllongreadsvimasm_vcf}.gz', f'{calllongreadsvimasm_vcfgz}')
                    shutil.rmtree(callgenomesvimasmtmp_path)
            except Exception as e:
                self.logger.error(f"An error occurred: {e}")

        return genomesvimasm_prefix, svimasm_log_out


    # Main method with multiprocessing
    def run_call(self):
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.jobs) as executor:
            futures = []
            call_path_dict = {}
            if self.fq_path:
                longreadminimap2_prefix = self.align_path_dict['longreadminimap2_prefix']
                calllongread_path = os.path.abspath(os.path.join(self.workdir_path, "02_call", "01_longread"))
                self.mkdirClass.make(calllongread_path, self.restart)

                future_longreadsniffles = executor.submit(self.run_sniffles_longread,
                                                          longreadminimap2_prefix, calllongread_path)
                futures.append(('longreadsniffles_prefix', future_longreadsniffles))
                future_longreadcuteSV = executor.submit(self.run_cuteSV_longread,
                                                          longreadminimap2_prefix, calllongread_path)
                futures.append(('longreadcuteSV_prefix', future_longreadcuteSV))
                future_longreaddebreak = executor.submit(self.run_debreak_longread,
                                                        longreadminimap2_prefix, calllongread_path)
                futures.append(('longreaddebreak_prefix', future_longreaddebreak))
                future_longreadsvim = executor.submit(self.run_svim_longread,
                                                         longreadminimap2_prefix, calllongread_path)
                futures.append(('longreadsvim_prefix', future_longreadsvim))
            if self.fa_path:
                genomeminimap2_prefix = self.align_path_dict['genomeminimap2_prefix']
                genomenucmer_prefix = self.align_path_dict['genomenucmer_prefix']
                callgenome_path = os.path.abspath(os.path.join(self.workdir_path, "02_call", "02_genome"))
                future_genomeAssemblytics = executor.submit(self.run_Assemblytics_genome,
                                                      genomenucmer_prefix, callgenome_path)
                futures.append(('genomeAssemblytics_prefix', future_genomeAssemblytics))
                future_genomesvimasm = executor.submit(self.run_svimasm_genome,
                                                            genomeminimap2_prefix, callgenome_path)
                futures.append(('genomesvimasm_prefix', future_genomesvimasm))

        for task_name, future in futures:
            call_path, log_out = future.result()
            if log_out:
                raise SystemExit(log_out)
            call_path_dict[task_name] = call_path

        return call_path_dict