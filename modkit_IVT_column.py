import pysam
import pickle
import pandas as pd
from tqdm import tqdm
from time import time
import argparse
from collections import defaultdict
    

def modkit_df(in_modkit_path):
    df = pd.read_csv(in_modkit_path, sep='\t', header=None)
    df.columns = ["chrom", "start", "stop", "mod", "score", "strand", "start2", "stop2", "color", "n_valid_cov", "proportion_modified",
                 "n_mod", "n_canon", "n_other_mod", "n_del", "n_fail", "n_diff", "n_nocall"]
    return df

def mod_to_col(mod_code_to_threshold, line):
    m_to_c = defaultdict(lambda: 2)
    for mod in mod_code_to_threshold:
        for i, l in enumerate(line.strip().split('\t')):
            if "false_positive" not in l:
                continue
            if 100 * float(l.split('_')[2]) > mod_code_to_threshold[mod]:
                if i == 2:
                    m_to_c[mod] = 2
                else:
                    m_to_c[mod] = i-1
                break
    return m_to_c

def load_error_table(error_table_path, mod_code_to_threshold):
    error_lookup_dict = {}
    with open(error_table_path, 'r') as infile:
        for i, line in enumerate(infile):
            if i == 0:
                mod_col = mod_to_col(mod_code_to_threshold, line)
            else:
                split_line = line.split('\t')
                kmer = split_line[0]
                mod = split_line[1]
                false_positive = split_line[mod_col[mod]]
                if mod not in error_lookup_dict:
                    error_lookup_dict[mod] = {}
                error_lookup_dict[mod][kmer] = float(false_positive)
    return error_lookup_dict

def ref_seq_dict(ref_path):
    seq_dict = {}
    ref_handle = pysam.FastaFile(ref_path)
    for ref in ref_handle.references:
        seq_dict[ref] = ref_handle.fetch(ref)
    return seq_dict

def main(in_modkit_path, 
         ref_path,
         error_table_path,
         mod_thresholds,
         outpath):

    mod_code_to_threshold = defaultdict(lambda: 0.7)
    for pair in mod_thresholds: 
        mod, threshold = pair.split(',')
        mod_code_to_threshold[mod] = float(threshold)
    
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = str.maketrans(old_chars,replace_chars)
    
    mod_to_base = {'17802':'T',
                   'a':'A',
                   'm':'C',
                   '17596':'A'}
    
    pre_df = modkit_df(in_modkit_path)
    error_table = load_error_table(error_table_path, mod_code_to_threshold)
    seq_dict = ref_seq_dict(ref_path)

#    for key in error_table:
#        for sub_key in error_table[key]:
#            print(f"{key} : {sub_key} : {error_table[key][sub_key]}")
    modulated_mod_occupancy = [0] * pre_df.shape[0]

    i = 0

    for chrom, pos, strand, mod, mod_occupancy in tqdm(zip(pre_df["chrom"].to_list(),
                                                           pre_df["start"].to_list(),
                                                           pre_df["strand"].to_list(),
                                                           pre_df["mod"].to_list(),
                                                           pre_df["proportion_modified"].to_list())):
        
        kmer = seq_dict[chrom][pos-4:pos+5]
        if strand == "-":
            kmer = kmer.translate(tab)[::-1]

        #print(f"{chrom=} {pos=} {mod=} {kmer=}")
        if len(kmer) != 9:
            modulated_mod_occupancy[i] = None
            i+= 1
            continue
        
        assert kmer[4] == mod_to_base[str(mod)]
        if kmer not in error_table[str(mod)]:
            modulated_mod_occupancy[i] = None
        else:
            modulated_mod_occupancy[i] = mod_occupancy - (100 * error_table[str(mod)][kmer])
        i += 1

    pre_df["mod_proportion-ivt"] = modulated_mod_occupancy
    print(f"{len(modulated_mod_occupancy)=}")
    print(f"{modulated_mod_occupancy[:10]=}")
    print(pre_df.head())
    pre_df.to_csv(outpath, sep="\t", header=None, index=False)
    
         

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--modkit",
                        "-m",
                        help="Modkit file to calculte the IVT controlled proportions for",
                        required=True)

    parser.add_argument("--reference",
                        "-ref",
                        help="Reference file used for alignment and modkit pileup",
                        required=True)

    parser.add_argument("--errortable",
                        "-e",
                        help="False positive rates for each 9mer",
                        required=True)

    parser.add_argument("--mod_threshold",
                        "-mt",
                        nargs="+",
                        help = "Mod code and threshold seperated by ',' for multiple modifications space separate the paired mod code and threshold",
                        required = False)

    parser.add_argument("--outpath",
                        "-o",
                        help="Outpath for IVT controlled modkit",
                        required=True)

    args = parser.parse_known_args()[0]

    main(args.modkit,
         args.reference,
         args.errortable,
         args.mod_threshold,
         args.outpath)
         
    
    #main(sys.argv[1], #modkit file
    #     sys.argv[2], #reference file
    #     sys.argv[3], #lookup_table
    #     sys.argv[4]) #New modkit path
         