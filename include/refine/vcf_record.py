#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"


import os
import gzip
import sys
import logging


# environment variables
env_path = {'PATH': os.environ.get('PATH')}
# include
include_dir = os.path.dirname(os.path.realpath(__file__))
for root, dirs, files in os.walk(os.path.join(include_dir, "include")):
    sys.path.append(root)

# VCF_record
class MyVcfRecord():
    # Flag whether global log handler has been added
    logger_initialized = False

    def __init__(self, vcflist):
        # If the log handler has not been added yet, add
        if not MyVcfRecord.logger_initialized:
            self.logger = logging.getLogger(__name__)
            formatter = logging.Formatter('[%(asctime)s] %(message)s')  # add date and time formatting strings
            handler = logging.StreamHandler()
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            MyVcfRecord.logger_initialized = True
        else:
            self.logger = logging.getLogger(__name__)

        # Value
        self.vcflist = vcflist
        self.vcfDict = {}

    # Build dict
    def build_dict(self):
        self.logger.error(f"Building vcf index: '{self.vcflist}' ...")
        for vcf_file in self.vcflist:
            with self.open_vcf_file(vcf_file) as f:
                for line in f:
                    if "#" in line:
                        continue
                    else:
                        line_list = line.strip().split("\t")
                        vcfId = line_list[2]
                        callerName = vcfId.split(".")[0]
                        self.vcfDict.setdefault(callerName, {})[vcfId] \
                            = line_list  # set the first string of the ID # column in the VCF file to the key of vcfDict

        return self.vcfDict

    # open file
    def open_vcf_file(self, vcf_file):
        if vcf_file.endswith('.gz'):
            return gzip.open(vcf_file, 'rt', encoding='UTF-8')
        else:
            return open(vcf_file, 'r', encoding='UTF-8')