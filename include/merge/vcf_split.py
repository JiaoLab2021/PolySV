#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"


import os
import sys
import re
import logging
import multiprocessing


# environment variables
env_path = {'PATH': os.environ.get('PATH')}
# include
include_dir = os.path.dirname(os.path.realpath(__file__))
for root, dirs, files in os.walk(os.path.join(include_dir, "include")):
    sys.path.append(root)

from vcf_io import MyVCF


class MyVcfsplit():
    # Flag whether global log handler has been added
    logger_initialized = False

    def __init__(self, threads, invcfpath, outvcfstandardpath, outvcflargepath, sequenceDict):
        # If the log handler has not been added yet, add
        if not MyVcfsplit.logger_initialized:
            self.logger = logging.getLogger(__name__)
            formatter = logging.Formatter('[%(asctime)s] %(message)s')  # add date and time formatting strings
            handler = logging.StreamHandler()  # output to the console
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            MyVcfsplit.logger_initialized = True
        else:
            self.logger = logging.getLogger(__name__)

        # Parameter
        self.threads = threads

        # Value
        self.invcfpath = invcfpath
        self.outvcfstandardpath = outvcfstandardpath
        self.outvcflargepath = outvcflargepath
        self.sequenceDict = sequenceDict

        # Class
        self.invcfClass = MyVCF(self.invcfpath)
        self.outvcfstandardClass = MyVCF(self.outvcfstandardpath)
        self.outvcflargeClass = MyVCF(self.outvcflargepath)


    # head
    def head(self):
        # Output
        output_head_lines = []
        invcfOpen = self.invcfClass.open_vcf_file()
        with invcfOpen as f_in:
            for line in f_in:
                if "#" in line:
                    output_head_lines.append(line)
                else:
                    break

        return output_head_lines


    # Split
    def split_sv(self, CHROMOSOME):

        # Storage List
        output_standard_lines = []
        output_large_lines = []

        # re
        patternTP = r"SVTYPE=(\w+)"
        patternSvlen = r"SVLEN=(-?\d+)"

        # I/O vcf
        invcfOpen = self.invcfClass.open_vcf_file()

        # Loop through
        with invcfOpen as f_in:
            for line in f_in:
                if "#" in line:
                    continue
                else:
                    info = line.strip()
                    infoList = info.split("\t")
                    INFOCol = infoList[7]
                    chromosome = infoList[0]

                    # Multi-process
                    if chromosome != CHROMOSOME:
                        continue

                    # Start and end position
                    refSeq = infoList[3]
                    svLen = abs(int(re.search(patternSvlen, INFOCol).group(1)))
                    svType = re.search(patternTP, INFOCol).group(1)

                    # DEL
                    if svType == "DEL":
                        if svLen < 100000:
                            output_standard_lines.append(line)
                        else:
                            output_large_lines.append(line)
                    elif svType == "INS":
                        if svLen < 20000:
                            output_standard_lines.append(line)
                        else:
                            output_large_lines.append(line)
                    elif svType == "INV":
                        svLen = len(refSeq) - 1
                        if svLen < 20000:
                            output_standard_lines.append(line)
                        else:
                            output_large_lines.append(line)
                    elif svType in ['DUP', 'Tandem_expansion', 'Tandem_contraction',
                                    'Repeat_expansion', 'Repeat_contraction']:
                        if svLen < 20000:
                            output_standard_lines.append(line)
                        else:
                            output_large_lines.append(line)

        return output_standard_lines, output_large_lines


    # result
    def output(self):
        # CHROMOSOME_list
        CHROMOSOME_list = list(self.sequenceDict.keys())  # Chromosome list

        # head result
        output_head_lines = self.head()
        # vcf result
        pool = multiprocessing.Pool(processes=self.threads)
        output_lines_pairs = pool.starmap(self.split_sv, [(CHROMOSOME,) for CHROMOSOME in CHROMOSOME_list])
        pool.close()
        pool.join()

        # out
        output_standard_lines = [pair[0] for pair in output_lines_pairs]
        output_large_lines = [pair[1] for pair in output_lines_pairs]
        outvcfstandardWrite = self.outvcfstandardClass.write_vcf_file()
        outvcflargeWrite = self.outvcflargeClass.write_vcf_file()
        with outvcfstandardWrite as f_standard:
            f_standard.writelines(output_head_lines)
            for output_standard_line in output_standard_lines:
                f_standard.writelines(output_standard_line)
        with outvcflargeWrite as f_large:
            f_large.writelines(output_head_lines)
            for output_large_line in output_large_lines:
                f_large.writelines(output_large_line)