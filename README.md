# LabUtilities
Collection of scripts that might be useful in daily lab work


######split_fasta.py
Splits a multi-sequence fasta file into multiple files containing n number of reads.
   Examples:
   Split a very large fasta file into smaller files containing 10,000 reads
     python split_fasta.py -f large_amino_acid_file.faa -n 10000
   NOTE: This script places the output files in the current directory.  It's recommended that
   you create an empty directory to work in.
