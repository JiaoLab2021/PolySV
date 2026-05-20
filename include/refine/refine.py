#!/usr/bin/env python3

# -*- coding: utf-8 -*-


__data__ = "2026/03/16"
__version__ = "1.0.6"
__author__ = "He Jiabao"
__email__ = "hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn"


import os
import re
import sys
import heapq
import logging
import multiprocessing
import math
import shutil
from datetime import datetime
from collections import Counter


# environment variables
env_path = {'PATH': os.environ.get('PATH')}
# include
include_dir = os.path.dirname(os.path.realpath(__file__))
for root, dirs, files in os.walk(os.path.join(include_dir, "include")):
    sys.path.append(root)

from make_dir import MyMkdir
from genome import MyGenome
from vcf_sort import MyVcfsort
from calculate_depth import MyDepth
from vcf_io import MyVCF


class MyRefine:
    # Flag whether global log handler has been added
    logger_initialized = False
    def __init__(self, merge_path, align_path_dict, vcfDict, sequenceDict, reference_path, fq_path, fa_path,
                 outvcf_path, samplename, jobs, threads, workdir_path, restart, l_depth_ratio, a_depth_ratio,
                 del_window, ins_window, command_line_str):
        # If the log handler has not been added yet, add
        if not MyRefine.logger_initialized:
            self.logger = logging.getLogger(__name__)
            formatter = logging.Formatter('[%(asctime)s] %(message)s')  # add date and time formatting strings
            handler = logging.StreamHandler()  # output to the console
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            MyRefine.logger_initialized = True
        else:
            self.logger = logging.getLogger(__name__)

        # Value
        self.merge_path = merge_path
        self.align_path_dict = align_path_dict
        self.vcfDict = vcfDict
        self.sequenceDict = sequenceDict

        # Parameter
        self.outvcf_path = outvcf_path
        self.reference_path = reference_path
        self.fq_path = fq_path
        self.fa_path = fa_path
        self.samplename = samplename
        self.workdir_path = workdir_path
        self.restart = restart
        # Algorithm
        self.l_depth_ratio = l_depth_ratio
        self.a_depth_ratio = a_depth_ratio
        self.del_window = del_window
        self.ins_window = ins_window
        self.inv_window = self.ins_window
        self.dup_window = self.del_window
        self.jobs = jobs
        self.threads = threads

        # Class
        self.referenceClass = MyGenome(self.reference_path)
        self.mergeClass = MyVCF(self.merge_path)
        self.outvcfClass = MyVCF(self.outvcf_path)
        self.mkdirClass = MyMkdir()

        # Variable
        self.longread_depth_dict = {}
        self.genome_depth_dict = {}
        self.longread_average_depth = ''
        self.genome_average_depth = ''

        # command_line
        self.command_line_str = command_line_str


    # Depth start breakpoint
    def depth_start_refine_part(self, start, depth_list, svCaller, longread_caller, window):
        # Calculate the depth of window above and below the starting breakpoint
        start_up = start - window
        start_down = start + window
        start_depth_list = depth_list[(start_up - 1):(start_down)]  # storage depth list
        # Correction of abnormal depth sites based on longread
        if svCaller in longread_caller:
            depth_ratio = self.l_depth_ratio
        # Correction of abnormal depth sites based on genome
        else:
            depth_ratio = self.a_depth_ratio
        # Update depth list and correction of abnormal depth sites based on genome
        start_depth_list[:3] = [sum(start_depth_list[:3]) / 3] * 3  # Modify the depth of the front of the list
        start_depth_list[-3:] = [sum(start_depth_list[-3:]) / 3] * 3  # Modify the depth of the end of the list
        for i in range(len(start_depth_list) - 4, 2, -1):
            if start_depth_list[i] < depth_ratio * (sum(start_depth_list[i + 1:i + 4]) / 3) \
                    and start_depth_list[i] < depth_ratio * (sum(start_depth_list[i - 3:i]) / 3):
                start_depth_list[i] = sum(start_depth_list[i:i + 3]) / 3
        # Whether the base is 0
        max_difference = None
        max_ratio = None
        if 0 in start_depth_list:
            dif_list = [start_depth_list[i] - start_depth_list[i + 1]
                        for i in range(len(start_depth_list) - 1)]  # Store depth difference list
            max_difference = float(max(dif_list))
            max_index = len(dif_list) - dif_list[::-1].index(max_difference) - 1  # Returns the last maximum index
        else:
            dif_list = [start_depth_list[i] / start_depth_list[i + 1]
                        for i in range(len(start_depth_list) - 1)]  # store depth ratio list
            max_ratio = float(max(dif_list))
            max_index = len(dif_list) - dif_list[::-1].index(max_ratio) - 1  # Returns the last maximum index

        return start_up, start_depth_list, max_index, max_difference, max_ratio


    # Depth end breakpoint
    def depth_end_refine_part(self, end, depth_list, svCaller, longread_caller, window):
        # Calculate the depth of window above and below the ending breakpoint
        end_up = end - window
        end_down = end + window
        end_depth_list = depth_list[(end_up - 1):(end_down)]  # storage depth list
        # Correction of abnormal depth sites based on longread
        if svCaller in longread_caller:
            depth_ratio = self.l_depth_ratio
        # Correction of abnormal depth sites based on genome
        else:
            depth_ratio = self.a_depth_ratio
        # Update depth list and correction of abnormal depth sites based on genome
        end_depth_list[:3] = [sum(end_depth_list[:3]) / 3] * 3  # Modify the depth of the front of the list
        end_depth_list[-3:] = [sum(end_depth_list[-3:]) / 3] * 3  # Modify the depth of the end of the list
        for i in range(3, len(end_depth_list) - 3, 1):
            if end_depth_list[i] < depth_ratio * (sum(end_depth_list[i + 1:i + 4]) / 3) \
                    and end_depth_list[i] < depth_ratio * (sum(end_depth_list[i - 3:i]) / 3):
                end_depth_list[i] = sum(end_depth_list[i - 3:i]) / 3
        # Whether the base is 0
        min_difference = None
        min_ratio = None
        if 0 in end_depth_list:
            dif_list = [end_depth_list[i] - end_depth_list[i + 1]
                        for i in range(len(end_depth_list) - 1)]  # Store depth difference list
            min_difference = float(min(dif_list))
            min_index = dif_list.index(min_difference)  # Returns the first minimum index
        else:
            dif_list = [end_depth_list[i] / end_depth_list[i + 1]
                        for i in range(len(end_depth_list) - 1)]  # store depth ratio list
            min_ratio = float(min(dif_list))
            min_index = dif_list.index(min_ratio)  # Returns the first minimum index

        return end_up, end_depth_list, min_index, min_difference, min_ratio


    # Output results
    def process_output_line(self, chromosome, refine_start_pos, refine_ref_seq, refine_qry_seq, supp,
                            infoList, svType, refine_end_pos, idTuplist, longread_caller):
        # CHROM
        chrom_list = chromosome
        # POS
        pos_list = str(refine_start_pos)
        # ID
        id_list = f'polysv.{svType}.{chromosome}.'
        # REF
        ref_list = refine_ref_seq
        # ALT
        alt_list = refine_qry_seq
        # QUAL
        qual_list = infoList[5]
        # FILTER
        filter_list = infoList[6]
        # INFO
        new_svLen = len(refine_qry_seq) - len(refine_ref_seq)
        info_list = f"SN={supp};SVTYPE={svType};SVLEN={new_svLen};END={refine_end_pos}"
        # FORMAT
        format_list = infoList[8]
        # SAMPLE
        tem_gt = []
        for idTup in idTuplist:
            if idTup[0] in longread_caller:
                tem_gt.append(self.vcfDict[idTup[0]][idTup[1]][9])  # store the genotype of each vcf record
        if tem_gt:  # this mutation record is supported by longread_caller
            gt_counts = Counter(tem_gt)
            sample_list = max(gt_counts, key=gt_counts.get)
        else:
            sample_list = '1/1'

        # Save
        tem_list = [chrom_list, pos_list, id_list, ref_list, alt_list, qual_list, filter_list, info_list, format_list,
                    sample_list]

        return tem_list


    # Identify variant types
    def svtype_info(self, INFOCol, svCaller, vcfId, patternIdlist, info, refStart, refEnd):
        # Infolist
        infoList = info.split("\t")

        # The ID of each variant's call record is used to determine the genotype.
        idList = re.search(patternIdlist, INFOCol).group(1).split(',')
        idTuplist = [(id.split(".")[0], id) for id in idList]

        # Keep the most appropriate merged results
        svtypeTuplist = [(element, element.split('.')[1]) for element in idList]
        # Normalize Assemblytics variant types
        svtypeTuplist = [(element[0], 'INS') if element[1] == 'Insertion' else
                         (element[0], 'DEL') if element[1] == 'Deletion' else
                         (element[0], 'DUP') if element[1] in {'Tandem_expansion', 'Tandem_contraction',
                                                               'Repeat_expansion', 'Repeat_contraction'} else
                         element for element in svtypeTuplist]
        svTypelist = [element[1] for element in svtypeTuplist]
        svType_counts = Counter(svTypelist)
        svType = max(svType_counts, key=svType_counts.get)
        # Renew vcfId
        Id = next((element[0] for element in svtypeTuplist if element[1] == svType), None)
        if Id != vcfId:
            vcfId = Id
            svCaller = vcfId.split('.')[0]
            infoList = self.vcfDict[svCaller][vcfId]
            Start = int(infoList[1])
            End = int(re.search(r"END=(\d+)", infoList[7]).group(1))
            infoList[7] = INFOCol.replace(f"END={refEnd}", f"END={End}")
            INFOCol = infoList[7]
        else:
            End = refEnd
            Start = refStart

        return idList, idTuplist, svType, vcfId, svCaller, infoList, svTypelist, Start, End, INFOCol


    # All DEL start and end breakpoint records
    def del_pos_index(self, CHROMOSOME):
        # Initialize the dictionary storing DEL sites
        del_pos_chr_dic = {key: [0, 0] for key in self.sequenceDict.keys() if key == CHROMOSOME}
        # re
        patternE = r"END=(\d+)"
        patternIdlist = r"IDLIST=([^;]+)"
        # I/O vcf
        mergevcfOpen = self.mergeClass.open_vcf_file()
        with mergevcfOpen as f:

            for line_vcf in f:
                info = line_vcf.strip()
                if "#" in info:  # comment line
                    continue
                else:
                    infoList = info.split("\t")
                    if len(infoList) <= 7:
                        print("BAD_LINE:", info)
                        continue
                    INFOCol = infoList[7]
                    vcfId = infoList[2][2:]
                    svCaller = vcfId.split('.')[0]
                    chromosome = infoList[0]
                    refStart = int(infoList[1])
                    refEnd = int(re.search(patternE, INFOCol).group(1))
                    # Multi-process
                    if chromosome != CHROMOSOME:
                        continue

                    idList, idTuplist, svType, vcfId, svCaller, infoList, svTypelist, refStart, refEnd, INFOCol = \
                        self.svtype_info(INFOCol, svCaller, vcfId, patternIdlist, info, refStart, refEnd)

                    if svType == 'DEL':
                        del_pos_chr_dic.setdefault(chromosome, []).append(refStart - 25)
                        del_pos_chr_dic.setdefault(chromosome, []).append(refEnd + 25)

        return del_pos_chr_dic


    # Depth Correction in INV
    def inv_depth_refine(self, svCaller, longread_caller, del_pos_list, depth_list, pos, window):
        # Average Depth
        small_average_window_depth = sum(depth_list[pos - 1 - window:pos - 1]) / window
        big_average_window_depth = sum(depth_list[pos:pos + window]) / window
        # Correction of abnormal depth sites based on longread
        if svCaller in longread_caller:
            depth_ratio = self.l_depth_ratio
        # Correction of abnormal depth sites based on genome
        else:
            depth_ratio = self.a_depth_ratio
        if depth_list[pos - 1] < depth_ratio * small_average_window_depth and \
                depth_list[pos - 1] < depth_ratio * big_average_window_depth:
            # Determine whether del_end is outside the range of del_pos_list
            if not any(i <= pos <= j for i, j in zip(del_pos_list[::2], del_pos_list[1::2])):
                depth_list[pos - 1] = (small_average_window_depth + big_average_window_depth) / 2

        return depth_list


    # Depth
    def depth_info(self):
        self.longread_depth_dict = {}
        self.genome_depth_dict = {}
        self.longread_average_depth = ''
        self.genome_average_depth = ''

        # path
        self.refine_tmppath = os.path.join(self.workdir_path, "05_refine", "tmp_work")
        self.mkdirClass.make(self.refine_tmppath, self.restart)

        self.lbam_path = ''
        self.abam_path = ''
        if self.fq_path:
            self.lbam_path = os.path.abspath(f'{self.align_path_dict["longreadminimap2_prefix"]}.bam')
        if self.fa_path:
            self.abam_path = os.path.abspath(f'{self.align_path_dict["genomeminimap2_prefix"]}.bam')
        # run
        self.depthClass = MyDepth(self.samplename, self.lbam_path, self.abam_path, self.refine_tmppath, self.jobs,
                                  self.threads, self.restart, self.sequenceDict)
        depth_path_dict = self.depthClass.run_depth()
        # Assignment
        if self.fq_path:
            self.longread_depth_dict = depth_path_dict['longreaddepth']
            longreaddepth_sum = 0
            longreaddepth_count = 0
            for depth in self.longread_depth_dict.values():
                longreaddepth_sum += sum(depth)
                longreaddepth_count += len(depth)
            self.longread_average_depth = longreaddepth_sum / longreaddepth_count if longreaddepth_count > 0 else 0
            self.logger.error(f'TGS alignment depth threshold: average depth {self.longread_average_depth}')
        if self.fa_path:
            self.genome_depth_dict = depth_path_dict['genomedepth']
            genomedepth_sum = 0
            genomedepth_count = 0
            for depth in self.genome_depth_dict.values():
                genomedepth_sum += sum(depth)
                genomedepth_count += len(depth)
            self.genome_average_depth = genomedepth_sum / genomedepth_count if genomedepth_count > 0 else 0
            self.genome_average_depth = math.ceil(self.genome_average_depth)
            self.logger.error(f'Genome alignment depth threshold: average depth {self.genome_average_depth}')


    # head
    def head(self):
        # initial
        initial_str = '##fileformat=VCFv4.3\n##source=polysv-1.0.6\n'
        # time
        current_time  = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        time_str = f'##fileDate={current_time }\n'
        # contig
        contigs = ''.join([f"##contig=<ID={name},length={len(leng)}>\n" for name, leng in self.sequenceDict.items()])
        # other
        other_str = '##ALT=<ID=DEL,Description="Deletion">\n' \
                    '##ALT=<ID=INS,Description="Insertion">\n' \
                    '##ALT=<ID=INV,Description="Inversion">\n' \
                    '##ALT=<ID=DUP,Description="Duplication">\n' \
                    '##FILTER=<ID=PASS,Description="All filters passed">\n' \
                    '##INFO=<ID=SN,Number=1,Type=Integer,Description="Number of samples supporting the variant">\n' \
                    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n' \
                    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n' \
                    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n' \
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        # command
        command_str = f'##conmmand={self.command_line_str}'
        # head
        head_str = f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.samplename}\n'
        output_head_lines = initial_str + time_str + contigs + other_str + command_str + head_str

        return output_head_lines


    # VCF refine
    def vcf_refine(self, CHROMOSOME):

        # Storage List
        output_lines = []

        # longread_caller
        longread_caller = ["sniffles", "svim", "cuteSV", "debreak"]
        # genome_caller
        genome_caller = ["svim-asm", "Assemblytics"]
        genome_caller_set = set(genome_caller)

        # re
        patternE = r"END=(\d+)"
        patternIdlist = r"IDLIST=([^;]+)"
        patternSupp = r"SUPP=(\d+)"

        # id
        nums_dic = {
            'delnum': 0,
            'invnum': 0,
            'dupnum': 0,
            'insnum': 0,
        }

        # Open vcf
        mergevcfOpen = self.mergeClass.open_vcf_file()

        # Loop through
        with mergevcfOpen as f:
            for line_vcf in f:
                info = line_vcf.strip()

                if "#" in info:
                    continue
                else:
                    infoList = info.split("\t")
                    INFOCol = infoList[7]
                    oldvcfId = infoList[2]
                    vcfId = infoList[2][2:]
                    oldsvType = vcfId.split(".")[1]
                    svCaller = vcfId.split('.')[0]
                    chromosome = infoList[0]
                    supp = re.search(patternSupp, INFOCol).group(1)
                    refStart = int(infoList[1])
                    refEnd = int(re.search(patternE, INFOCol).group(1))

                    # Multi-process
                    if chromosome != CHROMOSOME:
                        continue

                    idList, idTuplist, svType, vcfId, svCaller, infoList, svTypelist, refStart, refEnd, INFOCol = \
                        self.svtype_info(INFOCol, svCaller, vcfId, patternIdlist, info, refStart, refEnd)
                    refSeq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, refStart, refEnd)
                    refSeq_len = len(refSeq) - 1
                    qrySeq = infoList[4]

                    # Refined breakpoint
                    refine_start_pos = ''
                    refine_end_pos = ''
                    refine_ref_seq = ''
                    refine_qry_seq = ''

                    if svCaller in longread_caller:
                        depth_data = self.longread_depth_dict
                        average_depth = float(self.longread_average_depth)
                    else:
                        depth_data = self.genome_depth_dict
                        average_depth = float(self.genome_average_depth)
                    depth_list = depth_data[chromosome]
                    chromosome_length = len(depth_list)

                    # DEL

                    if svType == 'DEL':

                        #### start breakpoint

                        # start position
                        start = refStart
                        # Filter variant records with abnormal start and end sites
                        if start <= self.del_window:
                            self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:'
                                              f'{oldvcfId} -> Abnormal start site, remove ...')
                            continue
                        # Determine whether the depth at the original starting breakpoint is 0, and update the starting breakpoint value
                        del_filter = False
                        if (self.del_window + 1) <= start <= chromosome_length - self.del_window:
                            while start >= (self.del_window + 1) and depth_list[start - 1] < 0.5 * average_depth:
                                start -= 1
                                # Break out of inner loop
                                if start == self.del_window or refStart - start >= refSeq_len:
                                    del_filter = True
                                    break
                        else:
                            del_filter = True

                        # refine start breakpoint
                        start_up, start_depth_list, max_index, max_difference, max_ratio= self.depth_start_refine_part \
                            (start, depth_list, svCaller, longread_caller, self.del_window)

                        # Filter out abnormal mutation sites
                        if del_filter:
                            refine_start_pos = refStart
                        else:
                            refine_start_pos = start_up + max_index
                        if (max_ratio is not None and max_ratio < (100 / 80)) or \
                                (max_difference is not None and max_difference < (0.2 * average_depth)):
                            self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId}, '
                                              f'start breakpoint is no depth difference, remove ...')
                            continue

                        #### end breakpoint

                        # end position
                        end = refEnd
                        # Filter variant records with abnormal start and end sites
                        if end + self.del_window >= chromosome_length:
                            self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:'
                                              f'{oldvcfId} -> Abnormal end site, remove ...')
                        # Determine whether the depth at the original ending breakpoint is 0, and update the ending breakpoint value
                        del_filter = False
                        if (self.del_window + 1) <= end <= chromosome_length - self.del_window:
                            while end <= (chromosome_length - self.del_window) and depth_list[end - 1] < 0.5 * average_depth:
                                end += 1
                                # Break out of inner loop
                                if end == self.del_window or refEnd - end >= refSeq_len:
                                    del_filter = True
                                    break
                        else:
                            del_filter = True
                        # refine end breakpoint
                        end_up, end_depth_list, min_index, min_difference, min_ratio = self.depth_end_refine_part \
                            (end, depth_list, svCaller, longread_caller, self.del_window)

                        # Filter out abnormal mutation sites
                        if del_filter:
                            refine_end_pos = refEnd
                        else:
                            refine_end_pos = end_up + min_index
                        if (min_ratio is not None and min_ratio > 0.8) or \
                                (min_difference is not None and -min_difference < (0.2 * average_depth)):
                            self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId}, '
                                              f'end breakpoint is no depth difference, remove ...')
                            continue

                        ## Filter out false positive results
                        if refine_start_pos >= refine_end_pos:
                            self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId}, '
                                              f'false positive results, remove ...')
                            continue
                        del_depth_list = depth_list[(refine_start_pos - 1):refine_end_pos]
                        del_average_depth = sum(del_depth_list) / len(del_depth_list)
                        if del_average_depth > average_depth:
                            self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId}, '
                                              f'DEL average depth is {del_average_depth}, which is more than the '
                                              f'average depth of the whole genome {average_depth}, remove ...')
                            continue

                        # Corrected reference genome and query genome sequences
                        refine_ref_seq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, refine_start_pos, refine_end_pos)
                        refine_qry_seq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, refine_start_pos, refine_start_pos)
                        del_tem_list = self.process_output_line(chromosome, refine_start_pos, refine_ref_seq,
                                                                refine_qry_seq, supp, infoList, svType, refine_end_pos,
                                                                idTuplist, longread_caller)
                        nums_dic['delnum'] += 1
                        final_vcfId = f'{del_tem_list[2]}{nums_dic["delnum"]}'
                        del_tem_list[2] = final_vcfId  # final_vcfId
                        # log
                        log_message = (f'{chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId} -> '
                                       f'{chromosome}:{refine_start_pos}-{refine_end_pos}:{svType}:{final_vcfId}')
                        if refStart == refine_start_pos and refEnd == refine_end_pos:
                            self.logger.error(f'[skip] {log_message}')
                        else:
                            self.logger.error(f'[refine] {log_message}')
                        # save
                        del_tem_str = '\t'.join(del_tem_list).strip() + "\n"
                        output_lines.append(del_tem_str)


                    # INS

                    elif svType == 'INS':
                        '''
                        Filter false positive results (in only 2 mutation records that support finding SVs
                        based on longread and mutation records whose depth suddenly increases and are mistakenly
                        identified as INS)
                        '''
                        # Corrected reference genome and query genome sequences
                        refine_start_pos = refStart
                        # Filter variant records with abnormal start and end sites
                        if refStart <= self.ins_window or refStart + self.ins_window >= chromosome_length:
                            self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:'
                                              f'{oldvcfId} -> Abnormal start site, remove ...')
                            continue
                        refine_end_pos = refStart
                        refine_ref_seq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, refine_start_pos,
                                                                      refStart)
                        refine_qry_seq = refine_ref_seq + qrySeq[1:]

                        if len(idTuplist) == 2:
                            # Set ins conditional discriminator
                            ins_filter = False
                            longread2_set = set(item[0] for item in idTuplist)  # caller sets
                            if not longread2_set.intersection(genome_caller_set):
                                # self.ins_window
                                ins_up = refine_start_pos - self.ins_window
                                ins_down = refine_start_pos + self.ins_window
                                ins_depth_list = depth_list[(ins_up - 1):(ins_down)]

                                if 0 in ins_depth_list:
                                    if sum(ins_depth_list) == 0:
                                        self.logger.error(
                                            f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId}, '
                                            f'depth is 0, remove ...')
                                        continue
                                    else:
                                        min_sec_depth = heapq.nsmallest(2, set(ins_depth_list))[
                                            -1]  # returns the 2nd smallest element
                                        ins_depth_list = [(min_sec_depth * 0.25) if depth == 0 else depth for depth in
                                                          ins_depth_list]  # re-zero the site with a depth of 0

                                for i in range(1, len(ins_depth_list) - 1):
                                    if ins_depth_list[i] <= 0.75 * ins_depth_list[i + 1] and \
                                            ins_depth_list[i] <= 0.75 * ins_depth_list[i - 1]:
                                        min_depth = min(ins_depth_list[i - 1], ins_depth_list[i + 1])
                                        ins_depth_list[i] = min_depth
                                # Start filtering
                                for i in range(len(ins_depth_list) - 1):
                                    if ins_depth_list[i] / ins_depth_list[i + 1] < 0.75:
                                        ins_filter = True
                                        break
                                if ins_filter:
                                    continue

                        ins_tem_list = self.process_output_line(chromosome, refine_start_pos, refine_ref_seq,
                                                                refine_qry_seq, supp,infoList, svType, refine_end_pos,
                                                                idTuplist, longread_caller)
                        nums_dic['insnum'] += 1
                        final_vcfId = f'{ins_tem_list[2]}{nums_dic["insnum"]}'
                        ins_tem_list[2] = final_vcfId  # final_vcfId
                        # log
                        log_message = (f'{chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId} -> '
                                       f'{chromosome}:{refine_start_pos}-{refine_end_pos}:{svType}:{final_vcfId}')
                        if refStart == refine_start_pos and refEnd == refine_end_pos:
                            self.logger.error(f'[skip] {log_message}')
                        else:
                            self.logger.error(f'[refine] {log_message}')
                        # save
                        ins_tem_str = '\t'.join(ins_tem_list).strip() + "\n"
                        output_lines.append(ins_tem_str)

                    # INV

                    elif svType == 'INV':

                        del_pos_list = self.del_pos_dic[chromosome]

                        #### start breakpoint
                        start = refStart
                        # Filter variant records with abnormal start and end sites
                        if start <= self.inv_window:
                            self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:'
                                              f'{oldvcfId} -> Abnormal start site, remove ...')
                            continue

                        ## INV start breakpoint or DEL start breakpoint

                        # DEL start breakpoint

                        if sum(depth_list[(start - self.inv_window):(start)]) >= \
                                sum(depth_list[(start + 1):(start + self.inv_window + 1)]):

                            # Refine DEL start breakpoint

                            start_up, start_depth_list, max_index, max_difference, max_ratio = self. \
                                depth_start_refine_part(start, depth_list, svCaller, longread_caller, self.inv_window)
                            del_start = start_up + max_index
                            refine_start_pos = del_start

                            # Find INV starting breakpoint

                            # Initialize INV starting breakpoint
                            inv_start = del_start + 1
                            # Start iterating
                            while True:
                                # Check if the end of the list is reached
                                if inv_start >= refEnd:
                                    break
                                # Check first condition
                                elif depth_list[inv_start - 1] >= (sum(depth_list[inv_start:inv_start + 3]) / 3):
                                    # Update INV start breakpoint
                                    inv_start += 1
                                else:
                                    break
                                # Correct depth
                                depth_list = self.inv_depth_refine(svCaller, longread_caller, del_pos_list,
                                                                   depth_list, inv_start, self.inv_window)
                            # refine Start
                            inv_start_up, inv_start_depth_list, min_index, min_difference, min_ratio = self. \
                                depth_end_refine_part(inv_start, depth_list, svCaller, longread_caller, self.inv_window)
                            inv_start = inv_start_up + min_index + 1

                        # INV start breakpoint

                        else:

                            # Refine INV start breakpoint

                            start_up, start_depth_list, min_index, min_difference, min_ratio = self. \
                                depth_end_refine_part(start, depth_list, svCaller, longread_caller, self.inv_window)
                            inv_start = start_up + min_index + 1
                            refine_start_pos = inv_start

                            # Find DEL starting breakpoint

                            # Initialize DEL starting breakpoint
                            del_start = inv_start - 1
                            # Start iterating
                            while True:
                                # Check if the start of the list is reached
                                if del_start - self.inv_window <= 1 or inv_start - del_start >= refSeq_len:
                                    break
                                elif depth_list[del_start - 1] >= (sum(depth_list[del_start - 4:del_start - 1]) / 3):
                                    del_start -= 1
                                else:
                                    break
                                # Correct depth
                                depth_list = self.inv_depth_refine(svCaller, longread_caller, del_pos_list,
                                                                   depth_list, del_start, self.inv_window)
                            # refine Start
                            del_start_up, del_start_depth_list, max_index, max_difference, max_ratio = self. \
                                depth_start_refine_part(del_start, depth_list, svCaller, longread_caller, self.inv_window)
                            del_start = del_start_up + max_index

                        #### end breakpoint
                        end = refEnd
                        # Filter variant records with abnormal start and end sites
                        if end + self.inv_window >= chromosome_length:
                            self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:'
                                              f'{oldvcfId} -> Abnormal end site, remove ...')
                            continue
                        ## INV start breakpoint or DEL start breakpoint

                        # INV end breakpoint

                        if sum(depth_list[(end - self.inv_window):(end)]) >= \
                                sum(depth_list[(end + 1):(end + self.inv_window + 1)]):

                            # Refine INV end breakpoint

                            end_up, end_depth_list, max_index, max_difference, max_ratio = self. \
                                depth_start_refine_part(end, depth_list, svCaller, longread_caller, self.inv_window)
                            inv_end = end_up + max_index
                            refine_end_pos = inv_end

                            # Find DEL ending breakpoint

                            # Initialize DEL starting breakpoint
                            del_end = inv_end + 1
                            # Start iterating
                            while True:
                                # Check if the end of the list is reached
                                if del_end + self.inv_window <= chromosome_length or del_end - inv_end >= refSeq_len:
                                    break
                                # Check first condition
                                elif depth_list[del_end - 1] >= (sum(depth_list[del_end:del_end + 3]) / 3):
                                    # Update DEL end breakpoint
                                    del_end += 1
                                else:
                                    break
                                # Correct depth

                                depth_list = self.inv_depth_refine(svCaller, longread_caller, del_pos_list,
                                                                   depth_list, del_end, self.inv_window)
                            # refine End
                            del_end_up, del_end_depth_list, min_index, min_difference, min_ratio = self. \
                                depth_end_refine_part(del_end, depth_list, svCaller, longread_caller,
                                                        self.inv_window)
                            del_end = del_end_up + min_index

                        # DEL end breakpoint

                        else:

                            # Refine DEL end breakpoint

                            end_up, end_depth_list, min_index, min_difference, min_ratio = self. \
                                depth_end_refine_part(end, depth_list, svCaller, longread_caller, self.inv_window)
                            del_end = end_up + min_index
                            refine_end_pos = del_end

                            # Find INV end breakpoint

                            # Initialize INV starting breakpoint
                            inv_end = del_end - 1
                            # Start iterating
                            while True:
                                # Check if the start of the list is reached
                                if inv_end >= refStart:
                                    break
                                # Check first condition
                                elif depth_list[inv_end - 1] >= (sum(depth_list[inv_end:inv_end + 3]) / 3):
                                    inv_end -= 1
                                else:
                                    break
                                # Correct depth
                                depth_list = self.inv_depth_refine(svCaller, longread_caller, del_pos_list,
                                                                   depth_list, inv_end, self.inv_window)
                            # refine End
                            inv_end_up, inv_end_depth_list, max_index, max_difference, max_ratio = self. \
                                depth_start_refine_part(inv_end, depth_list, svCaller, longread_caller,
                                                        self.inv_window)
                            inv_end = inv_end_up + max_index

                        #### Keep mutation records

                        if del_start < inv_start < inv_end < del_end:

                            # INVDEL1
                            if inv_start - del_start - 1 >= 50:
                                refine_invdel1_start_pos = del_start
                                refine_invdel1_end_pos = inv_start - 1

                                # Filter out false positive results
                                invdel1_depth_list = depth_list[(refine_invdel1_start_pos - 1):refine_invdel1_end_pos]

                                if not (sum(invdel1_depth_list) / len(invdel1_depth_list) > average_depth):
                                    refine_invdel1_ref_seq = self.referenceClass.seq_extr(self.sequenceDict, chromosome,
                                                                                          refine_invdel1_start_pos,
                                                                                          refine_invdel1_end_pos)
                                    refine_invdel1_qry_seq = self.referenceClass.seq_extr(self.sequenceDict, chromosome,
                                                                                          refine_invdel1_start_pos,
                                                                                          refine_invdel1_start_pos)
                                    invdel1_svLen = len(refine_invdel1_qry_seq) - len(refine_invdel1_ref_seq)
                                    if -invdel1_svLen < 20000:  # unreliable if larger than 20kb
                                        invdel1_vcfId = f'polysv.DEL.{chromosome}.'
                                        invdel1_info_list = f"SN={supp};SVTYPE=DEL;SVLEN={invdel1_svLen};" \
                                                            f"END={refine_invdel1_end_pos}"
                                        # save
                                        tem_list = [chromosome, str(refine_invdel1_start_pos), invdel1_vcfId,
                                                    refine_invdel1_ref_seq, refine_invdel1_qry_seq] + \
                                                   infoList[5:7] + [invdel1_info_list] + infoList[8:]
                                        nums_dic['delnum'] += 1
                                        final_vcfId = f'{tem_list[2]}{nums_dic["delnum"]}'  # vcfid
                                        tem_list[2] = final_vcfId
                                        # log
                                        self.logger.error(f'[refine] {chromosome}:{refStart}-{refEnd}:{oldsvType}:'
                                                          f'{oldvcfId} -> {chromosome}:{refine_invdel1_start_pos}-'
                                                          f'{refine_invdel1_end_pos}:INVDEL1:{final_vcfId}')
                                        tem_str = '\t'.join(tem_list).strip() + "\n"
                                        output_lines.append(tem_str)

                                # INV record
                                refine_start_pos = inv_start
                            else:
                                # INV record
                                refine_start_pos = del_start

                            # INVDEL2
                            if del_end - inv_end >= 50:
                                refine_invdel2_start_pos = inv_end
                                refine_invdel2_end_pos = del_end

                                # Filter out false positive results
                                invdel2_depth_list = depth_list[(refine_invdel2_start_pos - 1):refine_invdel2_end_pos]

                                if not (sum(invdel2_depth_list) / len(invdel2_depth_list) > average_depth):
                                    refine_invdel2_ref_seq = self.referenceClass.seq_extr(self.sequenceDict, chromosome,
                                                                                          refine_invdel2_start_pos,
                                                                                          refine_invdel2_end_pos)
                                    refine_invdel2_qry_seq = self.referenceClass.seq_extr(self.sequenceDict, chromosome,
                                                                                          refine_invdel2_start_pos,
                                                                                          refine_invdel2_start_pos)
                                    invdel2_svLen = len(refine_invdel2_qry_seq) - len(refine_invdel2_ref_seq)
                                    if -invdel2_svLen < 20000:  # unreliable if larger than 20kb
                                        invdel2_vcfId = f'polysv.DEL.{chromosome}.'
                                        invdel2_info_list = f"SN={supp};SVTYPE=DEL;SVLEN={invdel2_svLen};" \
                                                            f"END={refine_invdel2_end_pos}"
                                        tem_list = [chromosome, str(refine_invdel2_start_pos), invdel2_vcfId,
                                                    refine_invdel2_ref_seq, refine_invdel2_qry_seq] + \
                                                   infoList[5:7] + [invdel2_info_list] + infoList[8:]
                                        nums_dic['delnum'] += 1
                                        final_vcfId = f'{tem_list[2]}{nums_dic["delnum"]}'  # vcfid
                                        tem_list[2] = final_vcfId
                                        # log
                                        self.logger.error(f'[refine] {chromosome}:{refStart}-{refEnd}:{oldsvType}:'
                                                          f'{oldvcfId} -> {chromosome}:{refine_invdel2_start_pos}-'
                                                          f'{refine_invdel2_end_pos}:INVDEL2:{vcfId}')
                                        # save
                                        tem_str = '\t'.join(tem_list).strip() + "\n"
                                        output_lines.append(tem_str)

                                # INV record
                                refine_end_pos = inv_end

                            else:
                                # INV record
                                refine_end_pos = del_end

                        # INV

                        refine_ref_seq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, refine_start_pos, refine_end_pos)
                        # if refine_ref_seq:
                        inv_refine_ref_seq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, inv_start, inv_end)
                        if refine_start_pos == inv_start:
                            if refine_end_pos == inv_end:
                                refine_qry_seq = self.referenceClass.dna_reverse_complement(refine_ref_seq)
                            elif refine_end_pos == del_end:
                                refine_qry_seq = self.referenceClass.dna_reverse_complement(
                                    inv_refine_ref_seq)
                        elif refine_start_pos == del_start:
                            refine_qry_seq = refine_ref_seq[0] + self.referenceClass.dna_reverse_complement(
                                inv_refine_ref_seq)
                            refine_qry_seq = refine_ref_seq[0] + self.referenceClass.dna_reverse_complement(
                                inv_refine_ref_seq)


                        inv_tem_list = self.process_output_line(chromosome, refine_start_pos, refine_ref_seq,
                                                                refine_qry_seq, supp, infoList, svType, refine_end_pos,
                                                                idTuplist, longread_caller)
                        nums_dic['invnum'] += 1
                        final_vcfId = f'{inv_tem_list[2]}{nums_dic["invnum"]}'  # vcfid
                        inv_tem_list[2] = final_vcfId
                        # log
                        self.logger.error(f'[refine] {chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId} '
                                          f'-> {chromosome}:{refine_start_pos}-{refine_end_pos}'
                                          f':{svType}:{vcfId}')
                        # save
                        inv_tem_str = '\t'.join(inv_tem_list).strip() + "\n"
                        output_lines.append(inv_tem_str)

                    # DUP

                    elif svType == 'DUP':

                        # Modify the mutation type of Assembly
                        idTuplist_dup = [(id[0], id[1].replace('Tandem_expansion', 'DUP')
                                          .replace('Tandem_contraction', 'DUP')
                                          .replace('Repeat_expansion', 'DUP')
                                          .replace('Repeat_contraction', 'DUP')) for id in idTuplist]
                        idTuplist_dup = [(first, second) for first, second in idTuplist_dup if 'DUP' in second]
                        # sort idTuplist_dup
                        dup_order = ['svi-masm', 'Assemblytics', 'syri', 'mumandco', 'sniffles', 'cuteSV', 'svim', 'debreak']
                        idTuplist_dup = sorted(idTuplist_dup, key=lambda x: dup_order.index(x[0]))
                        svCaller = idTuplist_dup[0][0]
                        # select caller id
                        vcf_id_dup = [tpl[1] for tpl in idTuplist if tpl[0] == svCaller][0]
                        # Determine whether the call result needs to be changed
                        if vcfId != vcf_id_dup:
                            vcfId = vcf_id_dup
                            infoList = self.vcfDict[svCaller][vcfId]
                            Start = int(infoList[1])
                            End = int(re.search(r"END=(\d+)", infoList[7]).group(1))
                            infoList[7] = INFOCol.replace(f"END={refEnd}", f"END={End}")
                            refStart = Start
                            refEnd = End
                            # Modify vcfid in output results
                            vcfId.replace('Tandem_expansion', 'DUP') \
                                .replace('Tandem_contraction', 'DUP') \
                                .replace('Repeat_expansion', 'DUP') \
                                .replace('Repeat_contraction', 'DUP')


                            #### start breakpoint

                            # start position
                            start = refStart
                            if start <= self.dup_window:
                                self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:'
                                                  f'{oldvcfId} -> Abnormal start site, remove ...')
                                continue
                            # Determine whether the depth at the original starting breakpoint is 0, and update the starting breakpoint value
                            dup_filter = False
                            if (self.del_window + 1) <= start <= chromosome_length - self.del_window:
                                while start <= (chromosome_length - self.del_window) and depth_list[
                                    start - 1] < 0.5 * average_depth:
                                    start += 1
                                    # Break out of inner loop
                                    if start == self.del_window or refStart - start >= refSeq_len:
                                        dup_filter = True
                                        break
                            else:
                                dup_filter = True

                            # refine start breakpoint
                            start_up, start_depth_list, min_index, min_difference, min_ratio = self \
                                 .depth_end_refine_part(start, depth_list, svCaller, longread_caller, self.del_window)

                            if dup_filter:
                                refine_start_pos = refStart
                            else:
                                refine_start_pos = start_up + min_index + 1
                            # Filter out abnormal mutation sites
                            if (min_ratio is not None and min_ratio > 0.8) or \
                                    (min_difference is not None and -min_difference < (0.2 * average_depth)):
                                self.logger.error(
                                    f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId}, '
                                    f'start breakpoint is no depth difference, remove ...')
                                continue

                            #### end breakpoint

                            # end position
                            end = refEnd
                            # Filter variant records with abnormal start and end sites
                            if end + self.dup_window >= chromosome_length:
                                self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:'
                                                  f'{oldvcfId} -> Abnormal end site, remove ...')
                            # Determine whether the depth at the original starting breakpoint is 0, and update the starting breakpoint value
                            dup_filter = False
                            if (self.del_window + 1) <= end <= chromosome_length - self.del_window:
                                while start <= (chromosome_length - self.del_window) and depth_list[
                                    start - 1] < 0.5 * average_depth:
                                    start += 1
                                    # Break out of inner loop
                                    if start == self.del_window or refStart - start >= refSeq_len:
                                        dup_filter = True
                                        break
                            else:
                                dup_filter = True
                            end_up, end_depth_list, max_index, max_difference, max_ratio = self \
                                .depth_start_refine_part(end, depth_list, svCaller, longread_caller, self.del_window)

                            if dup_filter:
                                refine_end_pos = refEnd
                            else:
                                refine_end_pos = end_up + max_index
                                # Filter out abnormal mutation sites
                                if (max_ratio is not None and max_ratio < (100 / 80)) or \
                                        (max_difference is not None and max_difference < (0.2 * average_depth)):
                                    self.logger.error(
                                        f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId}, '
                                        f'end breakpoint is no depth difference, remove ...')
                                    continue

                        else:

                            # refine start and end position
                            refine_start_pos = refStart
                            refine_end_pos = refEnd

                        # Save

                        # Filter out false positive results
                        if refine_start_pos >= refine_end_pos:
                            self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId}, '
                                              f'false positive results, remove ...')
                            continue
                        dup_depth_list = depth_list[(refine_start_pos - 1):(refine_end_pos)]
                        dup_average_depth = sum(dup_depth_list) / len(dup_depth_list)
                        if dup_average_depth < average_depth:
                            self.logger.error(f'[warning] {chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId}, '
                                              f'DUP average depth is {dup_average_depth}, which is less than the '
                                              f'average depth of the whole genome {average_depth}, remove ...')
                            continue
                        # Corrected reference genome and query genome sequences
                        refine_ref_seq = self.referenceClass.seq_extr(self.sequenceDict, chromosome, refine_start_pos,
                                                                      refine_end_pos)
                        refine_qry_seq = refine_ref_seq + "N" * 28 + refine_ref_seq

                        dup_tem_list = self.process_output_line(chromosome, refine_start_pos, refine_ref_seq,
                                                                refine_qry_seq, supp,infoList, svType, refine_end_pos,
                                                                idTuplist, longread_caller)
                        nums_dic['dupnum'] += 1
                        final_vcfId = f'{dup_tem_list[2]}{nums_dic["dupnum"]}'  # vcfid
                        dup_tem_list[2] = final_vcfId
                        # log
                        log_message = (f'{chromosome}:{refStart}-{refEnd}:{oldsvType}:{oldvcfId} -> '
                                       f'{chromosome}:{refine_start_pos}-{refine_end_pos}:{svType}:{final_vcfId}')
                        if refStart == refine_start_pos and refEnd == refine_end_pos:
                            self.logger.error(f'[skip] {log_message}')
                        else:
                            self.logger.error(f'[refine] {log_message}')
                        # save
                        dup_tem_str = '\t'.join(dup_tem_list).strip() + "\n"
                        output_lines.append(dup_tem_str)

        return output_lines


    # filter
    def filter_resule(self, output_lines):
        filtered_output_lines = []

        for output_line in output_lines:
            filtered_output_line = []

            # None Value
            previous_line = None
            for line in output_line:
                line_list = line.split("\t")

                if previous_line is None or line_list[:2] != previous_line[:2]:
                    filtered_output_line.append(line)
                previous_line = line_list

            filtered_output_lines.append(filtered_output_line)

        return filtered_output_lines

    def run_refine(self):
        # CHROMOSOME_list
        CHROMOSOME_list = list(self.sequenceDict.keys())  # Chromosome list

        # del_pos_index
        pool = multiprocessing.Pool(processes=self.threads)
        del_pos_chr_dics = pool.starmap(self.del_pos_index, [(CHROMOSOME,) for CHROMOSOME in CHROMOSOME_list])
        pool.close()
        pool.join()
        self.del_pos_dic = {}
        for del_pos_chr_dic in del_pos_chr_dics:
            self.del_pos_dic.update(del_pos_chr_dic)

        # depth_info
        self.depth_info()

        # head result
        output_head_lines = self.head()

        # vcf_refine
        self.logger.error(f'Start breakpoint optimization ...')
        pool = multiprocessing.Pool(processes=self.threads)
        output_lines = pool.starmap(self.vcf_refine, [(CHROMOSOME,) for CHROMOSOME in CHROMOSOME_list])
        pool.close()
        pool.join()
        self.logger.error(f'Breakpoint refined is completed ...')

        # filter result
        self.logger.error(f'Remove redundant variants records ...')
        filtered_output_lines = self.filter_resule(output_lines)
        self.logger.error(f'Redundant variants removal is completed ...')

        # out
        outvcfWrite = self.outvcfClass.write_vcf_file()
        self.logger.error(f'Save resulte ...')
        with outvcfWrite as f_out:
            f_out.writelines(output_head_lines)
            for filtered_output_line in filtered_output_lines:
                f_out.writelines(filtered_output_line)

        # sort
        self.sortClass = MyVcfsort(self.outvcf_path)
        self.sortClass.sort()

        # remove
        if os.path.exists(self.outvcf_path):
            shutil.rmtree(os.path.join(self.workdir_path, "05_refine"))