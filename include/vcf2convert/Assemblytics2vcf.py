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
from vcf_sort import MyVcfsort


class MyAssemblytics2convert:
    def __init__(self, reference_path, fa_path, samplename, threads, invcfpath, outvcfpath, sequenceDict):
        # Parameter
        self.reference_path = reference_path
        self.fa_path = fa_path
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
        self.faClass = MyGenome(self.fa_path)
        self.faDict = self.faClass.build_index()
        self.sortClass = MyVcfsort(self.outvcfpath)


    # head
    def head(self):
        # Head  VCF annotation information
        title = '##fileformat=VCFv4.3\n' \
                '##source=assemblytics\n' \
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n' \
                '##INFO=<ID=ChrB,Number=1,Type=String,Description="Chromoosme ID on the non-reference genome">\n' \
                '##INFO=<ID=StartB,Number=1,Type=Integer,Description="Start position on non-reference genome">\n' \
                '##INFO=<ID=EndB,Number=1,Type=Integer,Description="End position on non-reference genome">\n' \
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # Current running time
        for name, leng in self.referenceClass.sequenceDict.items():
            title += "##contig=<ID={},length={}>\n".format(name, len(leng))
        command_str = f"##polysv vcf2convertCommand=Assemblytics "  # Command
        header_line = command_str + time_str + "\n"
        title += header_line
        title += f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.samplename}\n'
        # Output
        output_head_lines = []
        output_head_lines.append(title)

        return output_head_lines


    # convert
    def convert(self, CHROMOSOME):

        # Storage List
        output_lines = []

        # id
        nums_dic = {
            'Deletionnum': 0,
            'Insertionnum': 0,
            'Tandem_expansionnum': 0,
            'Tandem_contractionnum': 0,
            'Repeat_expansionnum': 0,
            'Repeat_contractionnum': 0,
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

                    # Start and end position
                    refStart = int(infoList[1]) - 1
                    if refStart < 50:
                        continue
                    refEnd = int(infoList[2]) - 1
                    queryCol = re.split(':', infoList[9])
                    qryChr = queryCol[0]
                    qryStart = int(queryCol[1].split('-')[0]) - 1
                    qryEnd = int(queryCol[1].split('-')[1]) - 1
                    svType = infoList[6]
                    refGap = int(infoList[7])
                    qryGap = int(infoList[8])

                    # SVLEN
                    svLen = int(infoList[4])
                    if svLen < 50:  # Filter out SNPs and InDels
                        continue

                    # DEL
                    if "Deletion" == svType:
                        nums_dic['Deletionnum'] += 1
                        vcfId = f"Assemblytics.DEL.{CHROMOSOME}.{nums_dic['Deletionnum']}"
                        if qryGap < 0:
                            refStart = refStart + qryGap
                        svType = 'DEL'
                    # INS
                    elif "Insertion" in svType:
                        nums_dic['Insertionnum'] += 1
                        vcfId = f"Assemblytics.INS.{CHROMOSOME}.{nums_dic['Insertionnum']}"
                        if refGap < 0:
                            qryStart = qryStart + refGap
                        svType = 'INS'
                    # Tandem_expansion
                    elif "Tandem_expansion" in svType:
                        nums_dic['Tandem_expansionnum'] += 1
                        vcfId = f"Assemblytics.Tandem_expansion.{CHROMOSOME}.{nums_dic['Tandem_expansionnum']}"
                        if qryGap >= 0:
                            qryStart = qryStart + refGap
                            qryEnd = qryEnd - refGap
                        elif qryGap < 0:
                            qryStart = qryStart + refGap - qryGap
                            qryEnd = qryEnd - refGap + qryGap
                    # Tandem_contraction
                    elif "Tandem_contraction" in svType:
                        nums_dic['Tandem_contractionnum'] += 1
                        vcfId = f"Assemblytics.Tandem_contraction.{CHROMOSOME}.{nums_dic['Tandem_contractionnum']}"
                        if refGap >= 0:
                            refStart = refStart + qryGap
                            refEnd = refEnd - qryGap
                        elif refGap < 0:
                            refStart = refStart + qryGap - refGap
                            refEnd = refEnd - qryGap + refGap
                    # Repeat_expansion
                    elif "Repeat_expansion" in svType:
                        nums_dic['Repeat_expansionnum'] += 1
                        vcfId = f"Assemblytics.Repeat_expansion.{CHROMOSOME}.{nums_dic['Repeat_expansionnum']}"
                    # Repeat_contraction
                    elif "Repeat_contraction" in svType:
                        nums_dic['Repeat_contractionnum'] += 1
                        vcfId = f"Assemblytics.Repeat_contraction.{CHROMOSOME}.{nums_dic['Repeat_contractionnum']}"

                    ## Out
                    refSeq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, refStart, refEnd)
                    qrySeq = self.faClass.seq_extr(self.faDict, qryChr, qryStart, qryEnd)
                    if len(refSeq) == 0 or len(qrySeq) == 0:  # Filter out larger variations
                        continue
                    while refSeq[0] != qrySeq[0] and refStart > 1 and qryStart > 1:
                        refStart = refStart - 1
                        refSeq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, refStart, refStart) + refSeq
                        qryStart = qryStart - 1
                        qrySeq = self.faClass.seq_extr(self.faDict, qryChr, qryStart, qryStart) + qrySeq

                    # Final_line_str
                    INFOCol = f"SVTYPE={svType};END={refEnd};ChrB={qryChr};StartB={qryStart};EndB={qryEnd};" \
                              f"SVLEN={len(qrySeq) - len(refSeq)}"
                    Assemblytics_line_list = [chromosome, str(refStart), vcfId, refSeq, qrySeq, '.',
                                              'PASS', INFOCol, 'GT', '1/1']
                    final_line_str = '\t'.join(Assemblytics_line_list).strip() + "\n"
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
        # Sort
        self.sortClass.sort()

