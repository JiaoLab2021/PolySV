#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"

import io
import gzip
import logging


# Genome
class MyGenome():
    # Flag whether global log handler has been added
    logger_initialized = False

    def __init__(self, genome):
        # If the log handler has not been added yet, add
        if not MyGenome.logger_initialized:
            self.logger = logging.getLogger(__name__)
            formatter = logging.Formatter('[%(asctime)s] %(message)s')  # add date and time formatting strings
            handler = logging.StreamHandler()  # output to the console
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            MyGenome.logger_initialized = True
        else:
            self.logger = logging.getLogger(__name__)

        # Value
        self.genome = genome
        self.sequenceDict = {}


    # build index
    def build_index(self):
        self.logger.error(f"Building genome index: '{self.genome}' ...")
        chromosome = ""
        sequenceIO = io.StringIO()

        with self.open_genome_file() as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Not empty
                    if len(sequenceIO.getvalue()) > 0:
                        if not chromosome:
                            self.logger.error(f"Chromosome is empty ...")
                            exit(1)
                        # assignment
                        self.sequenceDict[chromosome] = sequenceIO.getvalue().strip()
                        self.logger.error(f"Length of {chromosome}: {len(self.sequenceDict[chromosome])} ...")
                        # clear memory
                        sequenceIO.seek(0)
                        sequenceIO.truncate(0)
                    chromosome = line.split()[0][1:]
                else:
                    sequenceIO.write(line)

        # Last save
        if len(sequenceIO.getvalue()) > 0:
            if not chromosome:
                self.logger.error(f"Chromosome is empty ...")
                exit(1)
            # Assignment
            self.sequenceDict[chromosome] = sequenceIO.getvalue().strip()
            self.logger.error(f"Length of {chromosome}: {len(self.sequenceDict[chromosome])} ...")
            # Alear memory
            sequenceIO.seek(0)
            sequenceIO.truncate(0)

        # Close the StringIO objects
        sequenceIO.close()

        return self.sequenceDict


    # Open file
    def open_genome_file(self):
        if self.genome.endswith('.gz'):
            return gzip.open(self.genome, 'rt', encoding='UTF-8')
        else:
            return open(self.genome, 'r', encoding='UTF-8')


    # Transform reverse complementary sequences
    def dna_reverse_complement(self, sequence):
        # Reverse sequence
        sequence = sequence[::-1]
        # Transform table
        trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhdNn', 'TGCAtgcaYRKMyrkmBVDHbvdhNn')
        string = sequence.translate(trantab)  # transform
        return string


    # Sequence extraction
    def seq_extr(self, sequenceDict, chromosome, start, end):
        if int(start) < int(end):
            seq = sequenceDict[chromosome][int(start) - 1: int(end)]
        else:
            seq = sequenceDict[chromosome][int(end) - 1: int(start)]
        return seq