# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#    File Name:      compute_frequency
#    Description:    The script is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered
#                    on-target efficiency.
#    Author:         Dengwei Zhang
#    Date:           2019/12/3
#    E-mail:         scdzzdw@163.com
#-------------------------------------------------
'''

import argparse
import os
import subprocess
import logging
import regex
import yaml
import textwrap


def parase_yaml(yaml_file, pam_position):
    with open(yaml_file, 'r') as f:
        items = yaml.load(f.read(), Loader=yaml.FullLoader)
        pcr_seq = items["pcr_seq"]
        target_seq = items["guide_rna"]["target"]
        pam_seq = items["guide_rna"]["PAM"]
        # define the pam_position is the end position for 5'-end PAM, or the start position for 3'-end PAM.
        if pam_position == '5_end':
            guide_seq = pam_seq + target_seq
            start_index = pcr_seq.find(guide_seq)
            if start_index == -1:
                logging.error("The guide RNA and PCR sequence are not on the same strand, please check!")
                os.exit(0)
            else:
                pam_end_position = start_index + len(pam_seq)
                end_index = start_index + len(guide_seq)
        else:
            guide_seq = target_seq + pam_seq
            start_index = pcr_seq.find(guide_seq)
            if start_index == -1:
                logging.error("The guide RNA and PCR sequence are not on the same strand, please check!")
                os.exit(0)
            else:
                pam_end_position = start_index + len(target_seq)
                end_index = start_index + len(guide_seq)
        adapter = items["adapter"]
    return pcr_seq, guide_seq, start_index, end_index, pam_end_position, adapter

def parase_cigar(Cigar, start_position):
	cigar_letter = [x for x in regex.split('\d', Cigar) if x != '']
	cigar_numer = [x for x in regex.split('\D', Cigar) if x != '']
	end_position = int(start_position)
	for i, k in zip(cigar_letter, cigar_numer):
		if i == "M" or i == "D":
			end_position += int(k)
		else:
			pass
	# end_position is the end position on the reference genome for a particular read
	return end_position

def regexFromSequence(seq, indels=20, substitution=6):
    seq = seq.upper()
    pattern = seq
    pattern_standard = '(' + '?b:' + pattern + ')' + '{{s<={0}}}'.format(substitution)
    pattern_gap = '(' + '?b:' + pattern + ')' + '{{i<={0},d<={0},s<={1}}}'.format(indels, substitution)
    return pattern_standard, pattern_gap

def alignDeletion(targetsite, deletion_seq, indels=20, substitution=6):
	pattern_standard, pattern_gap = regexFromSequence(targetsite, indels, substitution)
	m = regex.search(pattern_gap, deletion_seq)
	seq = deletion_seq
	for x in m.fuzzy_changes[2]:
		seq = seq[:x+1] + '-' + seq[x+1:]
	return seq

def alignSequences(targetsite, window_sequence, yaml_file, pam_position, left_most_position, indels=20, substitution=6):
    pcr_seq, guide_seq, start_index, end_index, pam_end_position, adapter = parase_yaml(yaml_file, pam_position)

    window_sequence = window_sequence.upper()
    query_regex_standard, query_regex_gap = regexFromSequence(targetsite, indels, substitution)

    alignments_mm = ('standard', regex.search(query_regex_standard, window_sequence, regex.BESTMATCH))
    alignments_bulge = ('gapped', regex.search(query_regex_gap, window_sequence, regex.BESTMATCH))

    final_alignment_seq, final_alignment_start, final_alignment_end, mutation_type  = "", "", ".", "."
    mutation_number, mutation_position = {"substitution":0, "deletion":0, "insertion":0}, {"substitution":"*", "deletion":"*", "insertion":"*"}
    
    # To recode each mutaion position and the number of each type of mutaion
    substitution_position, deletion_position, insertion_position = [], [], []
    substitution_number, deletion_number, insertion_number = 0, 0, 0

    # mutation_type: S for substitution, D for deletion, I for insertion, C for combined mutation

    alignment_type_m, match_m = alignments_mm

    if match_m != None:
        substitution, insertion, deletion = match_m.fuzzy_counts
        # determine whether the substitution postion exceeds the vulnerable region
        for x in match_m.fuzzy_changes[0]:
            if start_index - left_most_position <= x <= end_index - left_most_position:
                substitution_position.append(x)
                substitution_number += 1
        if substitution_position:
            mutation_type = "S" + str(substitution_number)
            mutation_position = {"substitution":substitution_position, "deletion":"*", "insertion":"*"}
            mutation_number = {"substitution":substitution_number, "deletion":0, "insertion":0}
            final_alignment_seq = match_m.group()
            final_alignment_start = match_m.start()
            final_alignment_end = match_m.end()

    if not final_alignment_seq:
        alignment_type_b, match_b = alignments_bulge
        if match_b != None:
            substitution, insertion, deletion = match_b.fuzzy_counts
            for s in match_b.fuzzy_changes[0]:
                if start_index - left_most_position <= s <= end_index -  left_most_position:
                    substitution_position.append(s)
                    substitution_number += 1
            if insertion < 20:
                for i in match_b.fuzzy_changes[1]:
                    if start_index - left_most_position <= i <= end_index - left_most_position:
                        insertion_position.append(i)
                        insertion_number +=1
            if deletion < 20:
                for d in match_b.fuzzy_changes[2]:
                    if start_index - left_most_position <= d <= end_index - left_most_position:
                        deletion_position.append(d)
                        deletion_number += 1
            if substitution_number == 0 and deletion_number == 0 and insertion_number != 0:
                mutation_type = "I" + str(insertion_number)
            elif substitution_number == 0 and insertion_number == 0 and deletion_number != 0:
                mutation_type = "D" + str(deletion_number)
            else:
                mutation_type = "C" + str(substitution_number + deletion_number + insertion_number)

            mutation_position = {"substitution":substitution_position, "deletion":match_b.fuzzy_changes[2], "insertion":match_b.fuzzy_changes[1]}
            mutation_number = {"substitution":substitution_number, "deletion":deletion, "insertion":insertion}
            seq = match_b.group()
            final_alignment_seq = alignDeletion(targetsite, seq, indels=20, substitution=6) if deletion else seq
            final_alignment_start = match_b.start()
            final_alignment_end = match_b.end()

    return [final_alignment_seq, final_alignment_start, final_alignment_end, mutation_type, mutation_position, mutation_number]

def fastp(read_1, read_2, yaml_file, pam_position):
    pcr_seq, guide_seq, start_index, end_index, pam_end_position, adapter = parase_yaml(yaml_file, pam_position)
    with open("adapter.fa", 'w') as f:
        for k, j in adapter.items():
            print(">", k, file=f)
            print(j, file=f)
    command = "fastp -i " + read_1 + " -o unmerged_read_1.fq.gz -I " + read_2 + " -O unmerged_read_2.fq.gz --thread 6 --merge " \
            "--merged_out merged_reads.fq.gz --adapter_fasta adapter.fa > log.fastp 2>&1"
    subprocess.check_call(command, shell=True)

def bowtie2(out_path, yaml_file, pam_position):
    pcr_seq, guide_seq, start_index, end_index, pam_end_position, adapter = parase_yaml(yaml_file, pam_position)
    bowtie2_path = os.path.join(out_path, "02.bowtie2")
    index_path = os.path.join(bowtie2_path, "index")
    alignment_path = os.path.join(bowtie2_path, "alignment")
    os.mkdir(alignment_path)
    os.mkdir(index_path)
    os.chdir(index_path)
    with open("pcr_product.fa", 'w') as f:
        print(">PCR_seq", file=f)
        print(pcr_seq, file=f)
    command_index = 'bowtie2-build pcr_product.fa index > log.index 2>&1'
    subprocess.check_call(command_index, shell=True)
    os.chdir(alignment_path)
    command_alignment = "bowtie2 --very-sensitive-local --threads 8 -x ../index/index -U ../../01.fastp/merged_reads.fq.gz -S align.sam > log.bowtie2 2>&1"
    subprocess.check_call(command_alignment, shell=True)

def parase_sam(out_path, yaml_file, pam_position, indels=20, substitution=6):
    sam_file = os.path.join(out_path, "02.bowtie2/alignment/align.sam")
    pcr_seq, guide_seq, start_index, end_index, pam_end_position, adapter = parase_yaml(yaml_file, pam_position)
    
    # recorde position distribution
    substitution_pos_distribution, deletion_pos_distribution, insertion_pos_distribution = {}, {}, {}
    mutation_type_distribution = {}
    total_clean_reads, mapped_reads, mutated_reads = 0, 0, 0
    
    with open(sam_file, 'r') as fin, open("intermediate.report", 'w') as fo, open("substitution_position_distribution", 'w') as fs, open("insertion_position_distribution", 'w') as fi, open("deletion_position_distribution", 'w') as fd, open("mutation_type_distribution", 'w') as fm, open("Statistics", 'w') as fsa:
        # print to intermediate.report
        print("ID\tDNA_Sequence\tMutation_type\tSubstitution_position\tDeletion_position\tInsertion_position\tSubstitution_number\tDeletion_number\tInsertion_number", file=fo)
        
        for line in fin:
            if not regex.search('^@', line):
                total_clean_reads += 1
                items = line.rstrip("\n").split("\t")
                read_name = items[0]
                flag = items[1]
                chromosome = items[2]
                leftmost_pos = int(items[3]) - 1 #  1-based leftmost mapping POSition in SAM file
                map_quality = items[4]
                cigar = items[5]
                mate_chromosome = items[6]
                mate_read_leftmost_pos = items[7]
                insert_size = items[8]
                sequence = items[9]
                if chromosome != "*":
                    mapped_reads += 1
                    reference_end_position = parase_cigar(cigar, leftmost_pos)
                    refernce_seq = pcr_seq[leftmost_pos:reference_end_position]

                    [final_alignment_seq, final_alignment_start, final_alignment_end, mutation_type, mutation_position, mutation_number] = alignSequences(refernce_seq, sequence, yaml_file, pam_position, leftmost_pos, indels, substitution)

                    if final_alignment_seq:
                        # re-define the coordinates in order to making its position be relative to PAM
                        for k, j in mutation_position.items():
                            if j != "*" and j:
                                if k == "substitution":
                                    mutation_position[k] = ','.join([str(x + leftmost_pos - pam_end_position + 1) for x in j])
                                    for x in j:
                                        correct_position = str(x + leftmost_pos - pam_end_position + 1)
                                        substitution_pos_distribution[correct_position] = substitution_pos_distribution.get(correct_position, 0) + 1
                                elif k == "insertion":
                                    mutation_position[k] = ','.join([str(x + leftmost_pos - pam_end_position + 1) for x in j])
                                    for x in j:
                                        correct_position = str(x + leftmost_pos - pam_end_position + 1)
                                        insertion_pos_distribution[correct_position] = insertion_pos_distribution.get(correct_position, 0) + 1
                                elif k == "deletion":
                                    mutation_position[k] = ','.join([str(x + leftmost_pos - pam_end_position + 2) for x in j])
                                    for x in j:
                                        correct_position = str(x + leftmost_pos - pam_end_position + 1)
                                        deletion_pos_distribution[correct_position] = deletion_pos_distribution.get(correct_position, 0) + 1
                            elif not j:
                                mutation_position[k] = "*"
                        mutation_type_distribution[mutation_type] = mutation_type_distribution.get(mutation_type, 0) + 1
                        print(read_name, '\t', final_alignment_seq, '\t', mutation_type, '\t', end="", file=fo)
                        print(mutation_position["substitution"], '\t', mutation_position["deletion"], '\t', mutation_position["insertion"], '\t', end="", file=fo)
                        print(mutation_number["substitution"], '\t', mutation_number["deletion"], '\t', mutation_number["insertion"], file=fo)
                        mutated_reads += 1

        # print distribution
        print("Position\tCounts", file=fs)
        print("Position\tCounts", file=fi)
        print("Position\tCounts", file=fd)
        print("Mutation_type\tCounts", file=fm)
        for k in sorted(substitution_pos_distribution.keys()):
            print(k, '\t', substitution_pos_distribution[k], file=fs)
        for k in sorted(insertion_pos_distribution.keys()):
            print(k, '\t', insertion_pos_distribution[k], file=fi)
        for k in sorted(deletion_pos_distribution.keys()):
            print(k, '\t', deletion_pos_distribution[k], file=fd)
        for k, v in sorted(mutation_type_distribution.items(),key=lambda x:x[1], reverse=True):
            print(k, '\t', v, file=fm)
        print("Total_clean_reads:\t", total_clean_reads, file=fsa)
        print("mapped_reads:\t", mapped_reads, file=fsa)
        print("mutated_reads:\t", mutated_reads, file=fsa)
        print("Mutation_rate:\t", '%.2f%%' % ((mutated_reads/mapped_reads) * 100), file=fsa)

def main():
    description = '''
------------------------------------------------------------------------------------------------------------------------

The script, only supporting paired-end reads temporarily, is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered on-target efficiency.

This is an example for input file with YAML format. Please note that the pcr sequence and guide RNA sequence must be one same strand.
# test.yaml
pcr_seq: TCTGTCTGAAACGGTCCCTGGCTAAACTCCACCCATGGGTTGGCCAGCCTTGCCTTGACCAATAGCCTTGACAAGGCAAACTTGACCAATAGTCTTAGAGTATCCAGTGAGGCCAGGGGCC
guide_rna:
    target: GCCAGCCTTGCCTTGACCAATAG
    PAM: TGGGTTG
adapter:
    seq_1: GCCAGCCTTGCCTTGACCAATAG
    seq_2: GCCAGCCTTGCCTTGACCAATAG
    seq_3: GCCAGCCTTGCCTTGACCAATAG
    seq_4: GCCAGCCTTGCCTTGACCAATAG
    
------------------------------------------------------------------------------------------------------------------------
'''

    parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))

    parse.add_argument("--read_1", help="the FASTQ file of read 1", required=True)
    parse.add_argument("--read_2", help="the FASTQ file of read 2", required=True)
    parse.add_argument("--output_path", help="the output path", required=True)
    parse.add_argument("--yaml", help="the configured YAML file", required=True)
    parse.add_argument("--pam_position", help="the PAM position, 5' end or 3' end", choices=['5_end', '3_end'], required=True)
    parse.add_argument("--indels", help="the maximum number of indels, and the maximum: 25, default: 20", default=20, type=int)
    parse.add_argument("--substitution", help="the maximum number of substitution, default: 6", default=6, type=int)

    args = parse.parse_args()
    read_1_path = os.path.abspath(args.read_1)
    read_2_path = os.path.abspath(args.read_2)

    # define the format of log
    log_file_path = os.path.join(args.output_path, 'my.log')
    LOG_FORMAT = "[%(asctime)s][%(levelname)s][%(module)s] %(message)s"
    DATE_FORMAT = "%m/%d/%Y %H:%M:%S %p"
    logging.basicConfig(filename=log_file_path, level=logging.DEBUG, format=LOG_FORMAT, datefmt=DATE_FORMAT)

    output_path = os.path.abspath(args.output_path)
    yaml_path = os.path.abspath(args.yaml)

    # run fastp
    logging.info("Running fastp...")
    fastp_path = os.path.join(output_path, "01.fastp")
    if not os.path.exists(fastp_path):
        os.mkdir(fastp_path)
    os.chdir(fastp_path)
    if not os.path.exists(os.path.join(fastp_path, "merged_reads.fq.gz")):
        fastp(read_1_path, read_2_path, yaml_path, args.pam_position)
    else:
        logging.warning("The mergered FASTQ has existed, skip fastp...")

    # run bowtie2
    logging.info("Running bowtie2...")
    bowtie2_path = os.path.join(output_path, "02.bowtie2")
    alignment_path = os.path.join(bowtie2_path, "alignment")
    sam_path = os.path.join(alignment_path, "align.sam")
    if not os.path.exists(bowtie2_path):
        os.mkdir(bowtie2_path)
    os.chdir(bowtie2_path)
    if not os.path.exists(sam_path):
        bowtie2(output_path, yaml_path, args.pam_position)
    else:
        logging.warning("The SAM file has existed, skip bowtie2...")

    # run analysis
    logging.info("Parasing SAM...")
    analysis_path = os.path.join(output_path, "03.analysis")
    if not os.path.exists(analysis_path):
        os.mkdir(analysis_path)
    os.chdir(analysis_path)
    parase_sam(output_path, yaml_path, args.pam_position, args.indels, args.substitution)

    logging.info("Finish analyszing...")


if __name__ == "__main__":
    main()