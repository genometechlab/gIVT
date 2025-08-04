#!/usr/bin/env python3
"""
Adjust modifications percentages in bedMethyl files by subtracting kmer-specific false positive rates.

This script reads an error table containing false positive rates for specific kmers,
then adjusts the modification percentages in a bedMethyl file accordingly.
"""

import pysam
import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
from collections import defaultdict
import gzip


def modkit_df(in_modkit_path):
    """
    Load bedMethyl file into pandas DataFrame with proper column names.
    
    Args:
        in_modkit_path: Path to bedMethyl file
        
    Returns:
        DataFrame with labeled columns for bedMethyl format
    """
    df = pd.read_csv(in_modkit_path, sep='\t', header=None)
    df.columns = ["chrom", "start", "stop", "mod", "score", "strand", "start2", "stop2", "color", 
                  "n_valid_cov", "proportion_modified", "n_mod", "n_canon", "n_other_mod", 
                  "n_del", "n_fail", "n_diff", "n_nocall"]
    return df


def mod_to_col(mod_code_to_threshold, line):
    """
    Determine which column contains the false positive rate for each modification.
    
    Args:
        mod_code_to_threshold: Dictionary mapping modification codes to thresholds
        line: Header line from error table
        
    Returns:
        Dictionary mapping modification codes to column indices
    """
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
    """
    Load error table from file, automatically handling gzipped files.
    
    Args:
        error_table_path: Path to error table (can be gzipped)
        mod_code_to_threshold: Dictionary of modification codes to thresholds
        
    Returns:
        Nested dictionary: {mod_code: {kmer: false_positive_rate}}
    """
    error_lookup_dict = {}
    
    # Determine if file is gzipped and open accordingly
    if error_table_path.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'  # Text mode for gzip
    else:
        open_func = open
        mode = 'r'
    
    with open_func(error_table_path, mode) as infile:
        for i, line in enumerate(infile):
            if i == 0:
                mod_col = mod_to_col(mod_code_to_threshold, line)
            else:
                split_line = line.strip().split('\t')
                kmer = split_line[0]
                mod = split_line[1]
                false_positive = split_line[mod_col[mod]]
                if mod not in error_lookup_dict:
                    error_lookup_dict[mod] = {}
                error_lookup_dict[mod][kmer] = float(false_positive)
    
    return error_lookup_dict


def ref_seq_dict(ref_path):
    """
    Load reference sequences into dictionary.
    
    This loads the entire reference into memory for fast random access.
    While this uses ~3.5GB for human genome, it's much faster than
    fetching sequences on-demand for millions of queries.
    
    Args:
        ref_path: Path to reference FASTA file
        
    Returns:
        Dictionary mapping chromosome names to sequences
    """
    seq_dict = {}
    ref_handle = pysam.FastaFile(ref_path)
    for ref in ref_handle.references:
        seq_dict[ref] = ref_handle.fetch(ref)
    ref_handle.close()  # Close handle after loading
    return seq_dict


def main(in_modkit_path, 
         ref_path,
         error_table_path,
         mod_thresholds,
         outpath):
    """
    Main function to adjust modifications percentages by subtracting false positive rates.
    
    Args:
        in_modkit_path: Path to input bedMethyl file
        ref_path: Path to reference genome
        error_table_path: Path to error table (can be gzipped)
        mod_thresholds: List of mod,threshold pairs
        outpath: Path for output bedMethyl file
    """
    # Parse modification thresholds
    mod_code_to_threshold = defaultdict(lambda: 0.7)
    if mod_thresholds:
        for pair in mod_thresholds: 
            mod, threshold = pair.split(',')
            mod_code_to_threshold[mod] = float(threshold)
    
    # Set up complement translation for reverse strand
    old_chars = "ACGT"
    replace_chars = "TGCA"
    tab = str.maketrans(old_chars, replace_chars)
    
    # Modification code to base mapping
    mod_to_base = {'17802': 'T',
                   '19229': 'G',
                   '19228': 'C',
                   '19227': 'T',
                   '69426': 'A',
                   'a': 'A',
                   'm': 'C',
                   '17596': 'A'}
    
    # Load data
    print("Loading bedMethyl file...")
    pre_df = modkit_df(in_modkit_path)
    
    print("Loading error table...")
    error_table = load_error_table(error_table_path, mod_code_to_threshold)
    
    print("Loading reference genome into memory...")
    seq_dict = ref_seq_dict(ref_path)
    
    # Get data as numpy arrays for better performance
    n_sites = len(pre_df)
    chroms = pre_df["chrom"].values
    positions = pre_df["start"].values
    strands = pre_df["strand"].values
    mods = pre_df["mod"].values
    occupancies = pre_df["proportion_modified"].values
    
    # Pre-allocate array for adjusted values
    adjusted_proportions = np.empty(n_sites, dtype=np.float64)

    # Process each modifications site
    print("Adjusting modifications percentages...")
    for i in tqdm(range(n_sites)):
        chrom = chroms[i]
        pos = positions[i]
        strand = strands[i]
        mod = str(mods[i])
        mod_occupancy = occupancies[i]
        
        # Extract 9mer context - using pre-loaded sequences for speed
        try:
            kmer = seq_dict[chrom][pos-4:pos+5]
        except (KeyError, IndexError):
            # Handle missing chromosomes or edge positions
            adjusted_proportions[i] = mod_occupancy
            continue
        
        # Reverse complement if on negative strand
        if strand == "-":
            kmer = kmer.translate(tab)[::-1]

        # Handle edge cases where kmer is too short
        if len(kmer) != 9:
            adjusted_proportions[i] = mod_occupancy
            continue
        
        # Handle case where central base doesn't match expected modification base
        if kmer[4] != mod_to_base.get(mod, 'N'):
            # Replace central base with modification code for lookup
            kmer = kmer[:4] + mod + kmer[5:]
            
        # Look up false positive rate and adjust
        if mod in error_table and kmer in error_table[mod]:
            # Subtract false positive rate and ensure non-negative
            adjusted_value = mod_occupancy - (100 * error_table[mod][kmer])
            adjusted_proportions[i] = max(0, adjusted_value)
        else:
            # Keep original value if kmer not in error table
            adjusted_proportions[i] = mod_occupancy

    # Replace original proportion_modified column with adjusted values
    pre_df["proportion_modified"] = adjusted_proportions
    
    # Write output in original bedMethyl format (no header, tab-separated)
    print("Writing output file...")
    pre_df.to_csv(outpath, sep="\t", header=None, index=False)
    
    # Print summary statistics
    print(f"\nProcessed {n_sites} modifications sites")
    print(f"Sample of adjusted values: {adjusted_proportions[:10]}")
    print(f"Output written to: {outpath}")
    
    # Get data as numpy arrays for better performance
    n_sites = len(pre_df)
    chroms = pre_df["chrom"].values
    positions = pre_df["start"].values
    strands = pre_df["strand"].values
    mods = pre_df["mod"].values
    occupancies = pre_df["proportion_modified"].values
    
    # Pre-allocate array for adjusted values
    adjusted_proportions = np.empty(n_sites, dtype=np.float64)

    # Process each modifications site
    print("Adjusting modifications percentages...")
    for i in tqdm(range(n_sites)):
        chrom = chroms[i]
        pos = positions[i]
        strand = strands[i]
        mod = str(mods[i])
        mod_occupancy = occupancies[i]
        
        # Extract 9mer context - using pre-loaded sequences for speed
        try:
            kmer = seq_dict[chrom][pos-4:pos+5]
        except (KeyError, IndexError):
            # Handle missing chromosomes or edge positions
            adjusted_proportions[i] = mod_occupancy
            continue
        
        # Reverse complement if on negative strand
        if strand == "-":
            kmer = kmer.translate(tab)[::-1]

        # Handle edge cases where kmer is too short
        if len(kmer) != 9:
            adjusted_proportions[i] = mod_occupancy
            continue
        
        # Handle case where central base doesn't match expected modification base
        if kmer[4] != mod_to_base.get(mod, 'N'):
            # Replace central base with modification code for lookup
            kmer = kmer[:4] + mod + kmer[5:]
            
        # Look up false positive rate and adjust
        if mod in error_table and kmer in error_table[mod]:
            # Subtract false positive rate and ensure non-negative
            adjusted_value = mod_occupancy - (100 * error_table[mod][kmer])
            adjusted_proportions[i] = max(0, adjusted_value)
        else:
            # Keep original value if kmer not in error table
            adjusted_proportions[i] = mod_occupancy

    # Add a adjusted proportion_modified column with adjusted values
    pre_df["adjusted_proportion_modified"] = adjusted_proportions
    
    # Write output in original bedMethyl format (no header, tab-separated)
    print("Writing output file...")
    pre_df.to_csv(outpath, sep="\t", header=None, index=False)

    print("Complete.")
         

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Adjust modifications percentages in bedMethyl files by subtracting kmer-specific false positive rates"
    )

    parser.add_argument("--modkit", "-m",
                        help="Modkit bedMethyl file to calculate the IVT-controlled proportions for",
                        required=True)

    parser.add_argument("--reference", "-ref",
                        help="Reference file used for alignment and modkit pileup",
                        required=True)

    parser.add_argument("--errortable", "-e",
                        help="False positive rates for each 9mer (can be gzipped)",
                        required=True)

    parser.add_argument("--mod_threshold", "-mt",
                        nargs="+",
                        help="Mod code and threshold separated by ',' for multiple modifications (e.g., m,0.7 a,0.8)",
                        required=False)

    parser.add_argument("--outpath", "-o",
                        help="Output path for IVT-controlled modkit bedMethyl file",
                        required=True)

    args = parser.parse_args()

    main(args.modkit,
         args.reference,
         args.errortable,
         args.mod_threshold,
         args.outpath)