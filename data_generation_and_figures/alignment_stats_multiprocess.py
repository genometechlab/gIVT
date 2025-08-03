import os
import sys
import pysam
import argparse
import numpy as np
from multiprocessing import Pool
import subprocess
from collections import defaultdict

def parse_cigar(cig):
    '''
    Return a list of cigar elements as tuples
    '''
    
    cigar_components = []
    current_num = ""
    for ch in cig:
        
        #split, add, reset on identified alpha char (I, D, M, H, S)
        if ch.isalpha() or ch == "=":
            if ch == "=" or ch == 'X':
                ch = 'M'
            cigar_components.append((int(current_num), ch))
            current_num = ""
        else:
            current_num += ch
    return cigar_components

def sum_cigar_components(parsed_cigar):
    n_count = 0; m_count = 0; s_count = 0; h_count = 0; d_count = 0; i_count = 0
    for element in parsed_cigar:
        if element[1] == 'M':
            m_count += element[0]
        elif element[1] == 'S':
            s_count += element[0]
        elif element[1] == 'H':
            h_count += element[0]
        elif element[1] == 'D':
            d_count += element[0]
        elif element[1] == 'I':
            i_count += element[0]
        elif element[1] == 'N':
            n_count += element[0]
    return (n_count, m_count, s_count, h_count, d_count, i_count)

def split_read_id_sets(bam_filepath, split_count):
    alignment_file = pysam.AlignmentFile(bam_filepath)
    list_of_sets = [set() for i in range(split_count)]
    
    i = 0
    count = 0
    for read in alignment_file.fetch():
        if (read.is_secondary or 
            read.is_unmapped or 
            read.is_supplementary):
            continue
        
        list_of_sets[i].add(read.query_name)
        count += 1
        i += 1
        if i >= len(list_of_sets):
            i = 0    
    
    return list_of_sets

def read_stats(pysam_read):
    
    cigar_stats = pysam_read.get_cigar_stats()
    n_count, m_count, s_count, h_count, d_count, i_count = sum_cigar_components(parse_cigar(pysam_read.cigarstring))
    
    
    read_id = pysam_read.query_name
    read_length = pysam_read.query_length
    contig = pysam_read.reference_name
    match = 0
    mismatch = 0
    
    for tup in pysam_read.get_aligned_pairs(with_seq=True):
        if tup[0] is None or tup[1] is None:
            continue
        if tup[2].islower():
            mismatch += 1
            continue
        match += 1
    align_length = match + d_count + mismatch
    read_identity = match / (i_count + d_count + match + mismatch)
    assert(match + mismatch == m_count)  
    return (read_identity, 
            read_id, 
            align_length, 
            match, 
            mismatch, 
            i_count, 
            d_count, 
            s_count, 
            h_count, 
            n_count,
            contig, 
            read_length)

def process_subset(args_list):
    bam_path = args_list[0]
    read_set = args_list[1]
    detailed_out_path = args_list[2]
    total_stats_out_path = args_list[3]
    
    read_identity_list = []
    align_length_list = []
    read_lengths = []
    
    running_total = {'aligned_bases':0,
                     'Insertions':0,
                     'Deletions':0,
                     'Matches':0,
                     'Mismatches':0,
                     'Soft Clipped':0,
                     'Skipped': 0,
                     'Hard Clipped':0
                    }
    
    infile = pysam.AlignmentFile(bam_path, 'rb')
    i = 0
    
    with open(detailed_out_path, 'w') as read_out_fh:
        for line in infile.fetch(until_eof=True):

            if (line.flag & 1 << 2 or 
                line.flag & 1 << 8 or 
                line.query_name not in read_set):
                continue


            stats = read_stats(line)
            read_lengths.append(stats[-1])
            running_total['aligned_bases'] += stats[2]
            running_total['Insertions'] += stats[5]
            running_total['Deletions'] += stats[6]
            running_total['Matches'] += stats[3]
            running_total['Mismatches'] += stats[4]
            running_total['Soft Clipped'] += stats[7]
            running_total['Skipped'] += stats[9]
            running_total['Hard Clipped'] += stats[8]

            read_identity_list.append(stats[0])
            align_length_list.append(stats[2])

            read_out_fh.write(f"{stats[1]}\t{stats[11]}\t{stats[10]}"+
                              f"\t{stats[2]}\t{stats[0]}\t{stats[3]}"+
                              f"\t{stats[4]}\t{stats[5]}\t{stats[6]}"+
                              f"\t{stats[7]}\t{stats[8]}\t{stats[9]}\n")
            i += 1
    
    with open(total_stats_out_path, 'w') as out_fh:
        read_identity_string = "\t".join([str(x) for x in read_identity_list])
        align_length_string = "\t".join([str(x) for x in align_length_list])
        read_lengths_string = "\t".join([str(x) for x in read_lengths])
                                         
        out_fh.write(f"{running_total['aligned_bases']}\n" +
                     f"{running_total['Insertions']}\n" +
                     f"{running_total['Deletions']}\n" +
                     f"{running_total['Matches']}\n" +
                     f"{running_total['Mismatches']}\n" +
                     f"{running_total['Soft Clipped']}\n" +
                     f"{running_total['Skipped']}\n" + 
                     f"{running_total['Hard Clipped']}\n" + 
                     f"{read_identity_string}\n" +
                     f"{align_length_string}\n" + 
                     f"{read_lengths_string}"
                    )
    return None

def main(args):
    
    read_id_sets = split_read_id_sets(args.input, args.threads)
    
    tmp_dir = os.path.join(args.output_directory, f"{args.output_prefix}_tmp_files")
    
    subprocess.run(["mkdir",
                    tmp_dir])
    
    args_list = [(args.input, 
                  read_id_sets[i], 
                  os.path.join(tmp_dir, f"tmp_read_stats_{i}.tsv"),
                  os.path.join(tmp_dir, f"tmp_overall_stats_{i}.tsv")) for i in range(len(read_id_sets))]
    
    with Pool(args.threads) as p:
        print(p.map(process_subset, args_list))
    
    running_total = defaultdict(int)
    
    read_stats_path = os.path.join(args.output_directory,
                                   f"{args.output_prefix}_read_stats.txt")

    gen_stats_path = os.path.join(args.output_directory,
                                  f"{args.output_prefix}_stats.txt")
    
    read_out_fh = open(read_stats_path, 'w')
    out_fh = open(gen_stats_path, 'w')
    read_out_fh.write("Read_ID\tRead_Length\tAligned_contig\t" + 
                      "Align_Length\tIdentity\tmatch\tmismatch" +
                      "\tinsertions\tdeletions\tsoftclip\thardclip" +
                      "\tintron_skip\n")
    
    read_identity_list = []
    align_length_list = []
    read_lengths_list = []
    
    for file in os.listdir(tmp_dir):
        fp = os.path.join(tmp_dir, file)
        
        if "read" in file:
            with open(fp, 'r') as infile:
                for line in infile:
                    read_out_fh.write(line)
        elif "tmp_overall_stats" in file:
            with open(fp, 'r') as infile:
                contents = infile.readlines()
                running_total["Matches"] += int(contents[3].strip())
                running_total["Mismatches"] += int(contents[4].strip())
                running_total["Deletions"] += int(contents[2].strip())
                running_total["Insertions"] += int(contents[1].strip())
                running_total["Skipped"] += int(contents[6].strip())
                running_total["aligned_bases"] += int(contents[0].strip())
                running_total["Hard Clipped"] += int(contents[7].strip())
                running_total["Soft Clipped"] += int(contents[5].strip())
                read_identity_list = read_identity_list + [float(x) for x in contents[8].strip().split("\t")]
                align_length_list = align_length_list + [int(x) for x in contents[9].strip().split("\t")]
                read_lengths_list = read_lengths_list + [int(x) for x in contents[10].strip().split("\t")]
        
    
    align_length_np = np.array(align_length_list)
    read_identity_np = np.array(read_identity_list)
    read_lengths_np = np.array(read_lengths_list)
    read_lengths_np.sort()
    read_lengths_np = np.flip(read_lengths_np)
    total_read_length = read_lengths_np.sum()
    cumsum_read_lengths = np.cumsum(read_lengths_np)
    print(f"{cumsum_read_lengths=}")
    print(f"{np.where(cumsum_read_lengths >= total_read_length / 2)[0]=}")
    print(f"{read_lengths_np[np.where(cumsum_read_lengths >= total_read_length / 2)[0]]=}")
    print(f"{read_lengths_np[np.where(cumsum_read_lengths >= total_read_length / 2)[0]][0]=}")
    n50 = read_lengths_np[np.where(cumsum_read_lengths >= total_read_length / 2)[0]][0]

    cumsum_aligned_read_lengths = np.cumsum(align_length_np)
    print(f"\n{align_length_np.sum()=}")
    print(f"{np.where(cumsum_read_lengths >= align_length_np.sum() / 2)=}")
    print(f"{np.where(cumsum_read_lengths >= align_length_np.sum() / 2)[0]=}")
    align_n50 = align_length_np[np.where(cumsum_aligned_read_lengths >= align_length_np.sum() / 2)[0]][0]
    print(running_total)
    print(align_length_np)
    print(read_identity_np)
    print(f"{n50=}")
    
    
    out_fh.write(f"Alignment Length Mean:                      {round(align_length_np.mean(), 4):18}\n")
    out_fh.write(f"Alignment Length Median:                    {round(np.median(align_length_np), 4):18}\n")
    out_fh.write(f"Read Length N50:                            {n50:18}\n")
    out_fh.write(f"Align Length N50:                           {align_n50:18}\n")

    out_fh.write(f"Read Identity Mean:                         {round(read_identity_np.mean(), 4):18}\n")
    out_fh.write(f"Read Identity Median:                       {round(np.median(read_identity_np), 4):18}\n")
    out_fh.write(f"Matches      per 1kb of reference bases:    {round(1000*running_total['Matches']/running_total['aligned_bases'],4):18}\n")
    out_fh.write(f"Mismatches   per 1kb of reference bases:    {round(1000*running_total['Mismatches']/running_total['aligned_bases'],4):18}\n")
    out_fh.write(f"Deletions    per 1kb of reference bases:    {round(1000*running_total['Deletions']/running_total['aligned_bases'],4):18}\n")
    out_fh.write(f"Insertions   per 1kb of reference bases:    {round(1000*running_total['Insertions']/running_total['aligned_bases'],4):18}\n")
    out_fh.write(f"Skipped      per 1kb of reference bases:    {round(1000*running_total['Skipped']/running_total['aligned_bases'],4):18}\n")


    out_fh.close()
    
    if args.detailed_read_stats:
        read_out_fh.close()
 
    
    subprocess.run(["rm",
                    "-rf",
                    tmp_dir])
                    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--input",
                        "-i",
                        required=True,
                        type=str,
                        help=("The input sam file for which alignment stats " +
                              "will be calculated.")
                       )
    
    parser.add_argument("--output_directory",
                        "-od",
                        required=True,
                        type=str,
                        help=("Path for directory that will contain output " +
                              "files")
                       )
        
    parser.add_argument("--output_prefix",
                        "-op",
                        required=True,
                        type=str,
                        help=("Prefix to append for output file identification")
                       )
    
    parser.add_argument("--overwrite",
                        default=False,
                        action="store_true")
    
    parser.add_argument("--threads",
                        "-t",
                        type = int,
                        default = 1)
    
    group = parser.add_mutually_exclusive_group()
    '''
    group.add_argument("--read_stats",
                        "-rs",
                        action="store_true",
                        help=("Records individual stats for each read, much " + 
                              "larger output file size.")
                       )
    '''
    group.add_argument("--detailed_read_stats",
                       "-drs",
                       action="store_true",
                       default=False,
                       help=("Records detailed individual read stats for " +
                             "each read. This includes the number of matches " +
                             "mismatches, insertions, deletions, softclip, " + 
                             "hardclip, and skipped (N, for Introns)")
                      )
    args = parser.parse_known_args()
    
    main(args[0])