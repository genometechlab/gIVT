# Modification False Positive Correction Tool

A Python tool for adjusting RNA modification measurements from Oxford Nanopore Dorado Basecaller by accounting for sequence-specific false positive rates.

## Overview

When measuring RNA modifications using nanopore sequencing, certain RNA sequences (kmers) can be incorrectly identified as modified even when they're not. This tool helps to correct for these systematic errors by subtracting sequence-specific false positive rates from your modification measurements.

## How It Works

The correction process follows these steps:

1. **Read your modification data**: The tool loads a bedMethyl file containing modification measurements for specific genomic positions, usually produced by modkit.
2. **Extract sequence context**: For each modified position, it looks up the surrounding 9-base DNA sequence (9mer) from the reference genome
3. **Look up error rates**: Using the provided error table, it finds the false positive rate for that specific 9mer sequence
4. **Apply correction**: It subtracts the false positive rate from the measured modification percentage
5. **Output corrected data**: The adjusted values are saved to a new bedMethyl file

## Installation
```bash
git clone https://github.com/genometechlab/gIVT
cd gIVT
pip install -r requirements.txt
chmod +x IVT_fp_correction.py
```

### Prerequisites

This tool requires Python 3.6 or higher and the following packages:

```bash
# Core requirements
pysam>=0.19.0
pandas>=1.3.0
numpy>=1.20.0
tqdm>=4.60.0
```

## Usage

### Basic Command
```bash
python IVT_fp_correction.py \
    --modkit input_bedMethyl.tsv \
    --reference reference_genome.fa \
    --errortable ./error_tables/GM12878_genomic_IVT_refmatch_9mer_all_threshold.dorado_1.0.tsv.gz
    --outpath corrected_bedMethyl.tsv
```

### Parameters

- **`--modkit` / `-m`**: Input bedMethyl file with modification data
- **`--reference` / `-ref`**: Reference genome FASTA (must match alignment reference)
- **`--errortable` / `-e`**: False positive rates per 9mer (can be gzipped)
- **`--outpath` / `-o`**: Output path for corrected bedMethyl file
- **`--mod_threshold` / `-mt`**: Optional modification thresholds (e.g., `m,0.7 a,0.8`), Default is 0.7
- **`--valid_only`**: Filters all output rows where a false positive adjustment couldn't be made
- **`--filter_mismatch`**: Filters all output rows where the original row doesn't match the reference base at the given position
- **`--filter_kmer`**: Filters all output rows where the representative kmer couldn't be found (too close to edge of contig, N in the reference).

### Example with Multiple Modifications

```bash
python apply_errors.py \
    --modkit multi_mod_methylation.tsv \
    --reference hg38.fasta \
    --errortable ./error_tables/GM12878_genomic_IVT_refmatch_9mer_all_threshold.dorado_1.0.tsv.gz \
    --mod_threshold m,0.7 a,0.8 \
    --outpath corrected_multimod.tsv \
    --valid_only
```

## Input File Formats
### bedMethyl
This program expects a bedMethyl format, specifically the same format produced by the modkit pileup program:
https://github.com/nanoporetech/modkit

The key columns are as follows:

- Column 1: Chromosome
- Column 2: Start position (0-based)
- Column 4: Modification code (e.g., 'm' for 5mC)
- Column 6: Strand (+/-)
- Column 11: Proportion modified (0-100)

Example:
```
chr1    10468    10469    m    0    +    10468    10469    255,0,0    33    6.06    2    31    0    0    0    0    0
chr1    10470    10471    m    0    +    10470    10471    255,0,0    33    9.09    3    30    0    0    0    0    0
```

### Error Table
The error table is a tsv file, each row is a specific 9mer and the columns represent modification and 9mer specific information:
```
mer             mod     false_positive_0.7_threshold    false_positive_0.71_threshold   ...   false_positive_0.99_threshold    Support
CTCATGGGC       17802   0.0100009185688697              0.009690901575345611                  0.0008037477609883801            87092 
```

Where:
- 'kmer': The 9-base sequence context (with base in question at the center)
- 'mod_code': The modification type output by Dorado Basecaller
- 'false_positive_{cutoff_threshold}': The observed rate of false positives given a specific mod quality cutoff threshold
- 'Support': Is the number of observed kmers in the original dataset that was used to calculate the false positive ratea

### Output 
The output file maintains the original structure, adding a false positive adjusted column as the 18th column (index:17)

### Notes


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

[Your acknowledgments here]
