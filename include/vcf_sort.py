#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"

# Import library
import os
import sys
import logging
import run_cmd


# environment variables
env_path = {'PATH': os.environ.get('PATH')}
# include
include_dir = os.path.dirname(os.path.realpath(__file__))
for root, dirs, files in os.walk(os.path.join(include_dir, "include")):
    sys.path.append(root)


# vcf
class MyVcfsort():
    # Flag whether global log handler has been added
    logger_initialized = False

    def __init__(self, vcfFile):
        # If the log handler has not been added yet, add
        if not MyVcfsort.logger_initialized:
            self.logger = logging.getLogger(__name__)
            formatter = logging.Formatter('[%(asctime)s] %(message)s')  # add date and time formatting strings
            handler = logging.StreamHandler()  # output to the console
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            MyVcfsort.logger_initialized = True
        else:
            self.logger = logging.getLogger(__name__)

        # Parameter
        self.vcfFile = vcfFile


    def sort(self):

        # Sort vcf file
        if self.vcfFile.endswith('.gz'):
            sortcmd = "zcat {0} | grep '#' | bgzip > {0}.sort; zcat {0} | grep -v '#' | " \
                "sort --parallel=10 -k 1,1 -k 2,2n -t $'\t' | bgzip >> {0}.sort; mv {0}.sort {0}".format(self.vcfFile)
        else:
            sortcmd = "grep '#' {0} > {0}.sort; grep -v '#' {0} | " \
                "sort --parallel=10 -k 1,1 -k 2,2n -t $'\t' >> {0}.sort; mv {0}.sort {0}".format(self.vcfFile)
        self.logger.error(f'Start sort the file {self.vcfFile}: {sortcmd}')
        run_cmd.run(sortcmd, env_path)
        self.logger.error(f'{self.vcfFile} sorting is done.')