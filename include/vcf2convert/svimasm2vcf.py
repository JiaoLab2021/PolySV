#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"


import os
import re
import sys
import multiprocessing
from datetime import datetime


# environment variables
env_path = {'PATH': os.environ.get('PATH')}
# include
include_dir = os.path.dirname(os.path.realpath(__file__))
for root, dirs, files in os.walk(os.path.join(include_dir, "include")):
    sys.path.append(root)

from genome import MyGenome
from vcf_io import MyVCF


class MySvimasm2convert:
    def __init__(self, reference_path, samplename, threads, invcfpath, outvcfpath, sequenceDict):

        # Parameter
        self.reference_path = reference_path
        self.samplename = samplename
        self.threads = threads

        # Value
        self.invcfpath = invcfpath
        self.outvcfpath = outvcfpath
        self.sequenceDict = sequenceDict

        # Class
        self.invcfClass = MyVCF(self.invcfpath)
        self.outvcfClass = MyVCF(self.outvcfpath)
        self.referenceClass = MyGenome(self.reference_path)


    # head
    def head(self):
        # Header Information
        time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # Current running time
        command_str = f"##polysv vcf2convertCommand=svim-asm "  # Command
        header_line = command_str + time_str + "\n"
        # Output
        output_head_lines = []
        invcfOpen = self.invcfClass.open_vcf_file()
        with invcfOpen as f_in:
            for line in f_in:
                info = line.strip()
                infoList = info.split("\t")
                if "##" in line:
                    output_head_lines.append(line)
                elif "#" in line:
                    output_head_lines.append(header_line)
                    line_list = infoList[:9] + [self.samplename]
                    line_str = '\t'.join(line_list).strip() + "\n"
                    output_head_lines.append(line_str)
                else:
                    break

        return output_head_lines


    # convert
    def convert(self, CHROMOSOME):
        # Storage List
        output_lines = []

        # re
        patternE = r"END=(\d+)"
        patternTP = r"SVTYPE=(\w+)"
        patternSvlen = r"SVLEN=(-?\d+)"

        # id
        nums_dic = {
            'delnum': 0,
            'invnum': 0,
            'dupnum': 0,
            'insnum': 0
        }

        # I/O vcf
        invcfOpen = self.invcfClass.open_vcf_file()

        # Loop through
        preline_info = None  # initialization
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

                    # Skip line
                    if "SVTYPE=BND" in infoList[7]:  # TRA
                        continue
                    if "SVTYPE=DUP" in infoList[7]:  # DUP
                        continue
                    if "PASS" not in infoList[6]:  # FILTER
                        continue

                    # Start and end position
                    refStart = int(infoList[1])
                    if refStart < 50:
                        continue
                    refEnd = int(re.search(patternE, INFOCol).group(1))
                    refSeq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, refStart, refEnd)
                    qrySeq = infoList[4]
                    svType = re.search(patternTP, INFOCol).group(1)

                    # SVLEN
                    if "SVTYPE=DEL" in INFOCol or "SVTYPE=INS" in INFOCol:
                        svLen = abs(len(qrySeq) - len(refSeq))
                    if "SVTYPE=INV" in INFOCol:
                        svLen = len(refSeq)
                    if svLen < 50:  # Filter out SNPs and InDels
                        continue

                    # Remove redundant SV
                    if preline_info:
                        preinfoList = preline_info.split("\t")
                        preINFOCol = preinfoList[7]
                        prechromosome = preinfoList[0]
                        prerefStart = int(preinfoList[1])
                        prerefEnd = int(re.search(patternE, preINFOCol).group(1))
                        presvType = re.search(patternTP, preINFOCol).group(1)
                        if (prechromosome == chromosome and
                                prerefStart == refStart and
                                prerefEnd == refEnd and
                                presvType == svType):
                            continue

                    # INV
                    if "SVTYPE=DEL" in INFOCol or "SVTYPE=INS" in INFOCol:
                        refSeq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, refStart, refEnd)
                        qrySeq = infoList[4]
                        SVLEN = int(re.search(patternSvlen, INFOCol).group(1))  # SVLEN value in line
                        INFOCol = INFOCol.replace(f"SVLEN={SVLEN}", f"SVLEN={len(qrySeq) - len(refSeq)}")
                        if "SVTYPE=DEL" in INFOCol:
                            nums_dic['delnum'] += 1
                            vcfId = f"svim-asm.DEL.{CHROMOSOME}.{nums_dic['delnum']}"
                        else:
                            nums_dic['insnum'] += 1
                            vcfId = f"svim-asm.INS.{CHROMOSOME}.{nums_dic['insnum']}"
                    elif "INV" == svType:
                        nums_dic['invnum'] += 1
                        refStart = refStart - 1
                        vcfId = f"svim-asm.INV.{CHROMOSOME}.{nums_dic['invnum']}"
                        refSeq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, refStart, refStart) + refSeq
                        qrySeq = refSeq[0] + qrySeq
                        INFOCol = INFOCol + ";SVLEN=" + str(len(qrySeq) - len(refSeq))

                    ## Out

                    # SAMPLE
                    genotype = re.split(':', infoList[9])[0]  # Genotype only (e.g. 1/1)

                    # Final_line_str
                    svimasm_line_list = [chromosome, str(refStart), vcfId, refSeq, qrySeq] + infoList[5:7] + \
                                        [INFOCol, "GT", genotype]
                    final_line_str = '\t'.join(svimasm_line_list).strip() + "\n"
                    output_lines.append(final_line_str)

                    # Update to a copy of the current row
                    preline_info = info


        return output_lines


    # result
    def output(self):
        # CHROMOSOME_list
        CHROMOSOME_list = list(self.sequenceDict.keys())  # Chromosome list

        # Head result
        output_head_lines = self.head()

        # VCF result
        pool = multiprocessing.Pool(processes=self.threads)
        output_lines = pool.starmap(self.convert, [(CHROMOSOME,) for CHROMOSOME in CHROMOSOME_list])
        pool.close()
        pool.join()

        # out
        outvcfWrite = self.outvcfClass.write_vcf_file()
        with outvcfWrite as f_out:
            f_out.writelines(output_head_lines)
            for output_line in output_lines:
                f_out.writelines(output_line)