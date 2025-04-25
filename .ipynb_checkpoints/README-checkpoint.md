# gIVT (Genomic IVT)
Whole Genome IVT for Nanopore DRS sequencing Control. This repository hosts a false positive error table for Oxford Nanopore Technologies modification caller Dorado. This error table is currently calculated for  Pseudouridine, m5C, m6A, and Inosine. Errors are calculated at a 9-mer level and are computed from [0.7-1.0) in increments of 0.01. A summary of the data and methods to produce this error table can be found [here](https://www.google.com).


## Controling for false positive calls
The modkit_IVT_column.py script can be run to add a column to the modkit output that subtracts the false positive rate for the associated modification and reference kmer. This tool does not work with modification calls at positions that do not match the reference kmer and may produce incorrect results at those positions. The minimum value of mod threshold calculated is 0.7, based on modkit recommendations this is the recommended. minimum confidence threshold

### Executing modkit_IVT_column.py
```python modkit_IVT_column.py --modkit {input_modkit_path} --reference {reference_path} --errortable {path_to_error_table} --outpath {output_path} --mod_threshold {mod_code,threshold} {mod_code,threshold}```
