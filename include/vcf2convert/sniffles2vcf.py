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


class MySniffles2convert:
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
        command_str = f"##polysv vcf2convertCommand=sniffles "  # Command
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
        patternR = r"SUPPORT=(\d+)"  # Read Depth
        patternE = r"END=(\d+)"
        patternL = r"SVLEN=(-?\d+)"

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
                    if "SVTYPE=BND" in info:  # TRA
                        continue
                    if int(infoList[5]) < 10:  # QUAL
                        continue
                    if infoList[6] != "PASS":  # FILTER
                        continue
                    if "IMPRECISE" in INFOCol:  # Imprecise breakpoint mutations
                        continue
                    readDepth = int(re.search(patternR, INFOCol).group(1))
                    if readDepth < 5:  # Read Depth
                        continue

                    # Start and end position
                    refStart = int(infoList[1]) - 1
                    if refStart < 50:
                        continue

                    refEnd = int(re.search(patternE, INFOCol).group(1)) - 1
                    #refSeq = self.referenceClass.seq_extr(self.sequenceDict, CHROMOSOME, refStart, refEnd)
                    qrySeq = infoList[4]

                    # SVLEN
                    svLen = abs(int(re.search(patternL, INFOCol).group(1)))
                    if "SVTYPE=INS" in INFOCol:
                        svLen = len(qrySeq)
                    if svLen < 50:  # Filter out SNPs and InDels
                        continue

                    # DEL
                    if "SVTYPE=DEL" in INFOCol:
                        nums_dic['delnum'] += 1
                        # 51bp to 100kb DEL
                        if svLen > 100000:
                            continue
                        vcfId = f"sniffles.DEL.{CHROMOSOME}.{nums_dic['delnum']}"
                        refSeq = self.referenceClass.seq_extr(self.sequenceDict, CHROMOSOME, refStart, refEnd)
                        qrySeq = self.referenceClass.seq_extr(self.sequenceDict, CHROMOSOME, refStart, refStart)
                    # INS
                    elif "SVTYPE=INS" in INFOCol:
                        nums_dic['insnum'] += 1
                        # 51bp to 20kb INS
                        if svLen > 20000 or "INS" in infoList[4]:
                            continue
                        vcfId = f"sniffles.INS.{CHROMOSOME}.{nums_dic['insnum']}"
                        refSeq = self.referenceClass.seq_extr(self.sequenceDict, CHROMOSOME, refStart, refStart)
                        qrySeq = refSeq + infoList[4]
                    # INV
                    elif "SVTYPE=INV" in INFOCol:
                        nums_dic['invnum'] += 1
                        # 51bp to 20kb INV
                        if svLen > 20000:
                            continue
                        vcfId = f"sniffles.INV.{CHROMOSOME}.{nums_dic['invnum']}"
                        refSeq = self.referenceClass.seq_extr(self.sequenceDict, CHROMOSOME, refStart, refEnd)
                        qrySeq = refSeq[0] + self.referenceClass.dna_reverse_complement(refSeq[1:])
                    # DUP
                    elif "SVTYPE=DUP" in INFOCol:
                        nums_dic['dupnum'] += 1
                        # 51bp to 20kb DUP
                        if svLen > 20000:
                            continue
                        vcfId = f"sniffles.DUP.{CHROMOSOME}.{nums_dic['dupnum']}"
                        refSeq = self.referenceClass.seq_extr(self.sequenceDict, CHROMOSOME, refStart, refEnd)
                        qrySeq = refSeq + "N"*28 + refSeq[1:]

                    ## OUT

                    # SVLEN
                    svLen = len(qrySeq) - len(refSeq)
                    INFOCol = re.sub(r"(?<=SVLEN=)-?\d+", str(svLen), INFOCol)

                    # END
                    INFOCol = re.sub(r"(?<=END=)\d+", str(refEnd), INFOCol)

                    # SAMPLE
                    genotype = re.split(':', infoList[9])[0]  # Genotype only (e.g. 1/1 or 0/1 or 1/0)
                    if genotype in ["./.", "0/0", ".", "0"]:
                        continue

                    # Final_line_str
                    sniffles_line_list = [CHROMOSOME, str(refStart), vcfId, refSeq, qrySeq] + infoList[5:7] + \
                              [INFOCol, "GT", genotype]
                    final_line_str = '\t'.join(sniffles_line_list).strip() + "\n"
                    output_lines.append(final_line_str)

        return output_lines


    # result
    def output(self):
        # CHROMOSOME_list
        CHROMOSOME_list = list(self.sequenceDict.keys())  # Chromosome list

        # head result
        output_head_lines = self.head()
        # vcf result
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