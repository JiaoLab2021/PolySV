#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"

import gzip

# VCF
class MyVCF():
    def __init__(self, vcfFile):
        # Vcf
        self.vcfFile = vcfFile


    # Open file
    def open_vcf_file(self):
        if self.vcfFile.endswith('.gz'):
            return gzip.open(self.vcfFile, 'rt', encoding='UTF-8')
        else:
            return open(self.vcfFile, 'r', encoding='UTF-8')

    # Write file
    def write_vcf_file(self):
        if self.vcfFile.endswith('.gz'):
            return gzip.open(self.vcfFile, 'wt', encoding='UTF-8')
        else:
            return open(self.vcfFile, 'w', encoding='UTF-8')