import pandas as pd
import pickle
import numpy as np
import pysam
import seaborn as sns
import sys
import time
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')

def mer_dict_to_df(mer_dict):
    total = 0
    for mer in mer_dict:
        for mod in mer_dict[mer]:
            total += 1
    false_positive = [0] * total
    mods = [0] * total
    mer_list = [0] * total
    support = [0] * total

    i = 0
    for mer in mer_dict:
        for mod in mer_dict[mer]:
            false_positive[i] = mer_dict[mer][mod]["mod"] / (mer_dict[mer][mod]["canon"] + mer_dict[mer][mod]["mod"])
            mods[i] = mod
            mer_list[i] = mer
            support[i] = mer_dict[mer][mod]["canon"] + mer_dict[mer][mod]["mod"]
            i += 1

    mer_dict_df = pd.DataFrame()
    mer_dict_df["false_positive"] = false_positive
    mer_dict_df["mod"] = mods
    mer_dict_df["mer"] = mer_list
    mer_dict_df["support"] = support
    mer_dict_df.sort_values(by="false_positive", ascending=False, inplace=True)
    mer_dict_df["x_axis"] = [i for i in range(len(false_positive))]

    

    return mer_dict_df

def reverse_compliment(intake, comp_dict): 
    return "".join([comp_dict[x] for x in list(reversed(intake))])

def read_id_lens(bam_path, min_len = 0, max_len = np.inf):
    af = pysam.AlignmentFile(bam_path)
    read_set = set()
    for read in af:
        if max_len >= read.query_alignment_length >= min_len:
            read_set.add(read.query_name)
    return read_set

def mer_mods(extract_full_path, mer_size, read_ids = None):
    mer_all = {}
    mer_ref_match = {}
    comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', '.':'.'}

    with open(extract_full_path, 'r') as infile:
        for line in tqdm(infile):
            if line.split('\t')[0] == "read_id":
                continue
            (read_id, 
            forward_read_position, 
            ref_position, 
            chrom, 
            mod_strand,
            ref_strand, 
            ref_mod_strand, 
            fw_soft_clipped_start, 
            fw_soft_clipped_end, 
            read_length, 
            mod_qual, 
            mod_code, 
            base_qual, 
            ref_kmer, 
            query_kmer, 
            cannonical_base, 
            modified_primary_base, 
            inferred, 
            flag) = line.strip().split('\t')  

            mod_qual = float(mod_qual)
            if read_ids is not None and read_id not in read_ids:
                continue
                
            if "-" in query_kmer:
                continue

            if sum([x in comp_dict for x in ref_kmer]) != mer_size:
                continue
            
            if ref_strand == '-':
                ref_kmer = reverse_compliment(ref_kmer, comp_dict)

            
            if mer_size < len(ref_kmer) or mer_size < len(query_kmer):
                chop_off = int((len(ref_kmer) - mer_size)/2)
                ref_kmer = ref_kmer[chop_off:-1 * chop_off]
                query_kmer = query_kmer[chop_off:-1 * chop_off]

            if mer_size < len(ref_kmer) or mer_size < len(query_kmer):
                print(f"{ref_kmer=}")
                x = input()
                if x == " ":
                    break
                    
            if ref_kmer == query_kmer:
                if query_kmer not in mer_ref_match:
                    mer_ref_match[query_kmer] = {}
                if mod_code not in mer_ref_match[query_kmer]:
                    mer_ref_match[query_kmer][mod_code] = {}
                    for i in range(30):
                        mer_ref_match[query_kmer][mod_code][round(0.7 + i / 100, 4)] = {'mod':0, 'canon':0}
                for fp_level in mer_ref_match[query_kmer][mod_code]:
                    if mod_qual >= fp_level:
                        mer_ref_match[query_kmer][mod_code][fp_level]["mod"] += 1
                    else:
                        mer_ref_match[query_kmer][mod_code][fp_level]["canon"] += 1
    return mer_ref_match

def main(inpath, outpath):
    extract_full_path = inpath

    ninemer_ref = mer_mods(extract_full_path, 9)
    with open(outpath, 'wb') as outfile:
        pickle.dump(ninemer_ref, outfile)
  
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])