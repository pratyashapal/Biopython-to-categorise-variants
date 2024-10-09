This Python script processes a VCF file to analyze genetic variants. 
It checks the quality of each variant, determines if it's in a coding region, and assesses its impact on the protein sequence. 
The script takes a VCF file, a GFF annotation file, and a genome FASTA file as inputs via command line.

Features

Variant Quality Check: Filters variants based on quality (QUAL > 20).

Coding Region Analysis: Identifies if a variant falls within a coding region.

Protein Impact: Determines if coding region variants are synonymous (no amino acid change) or non-synonymous (causes a change).

Error Logging: Logs errors and issues in processing files.

Outputs

1. Log File: Reports:

Input file names.

Count of low-quality variants (QUAL ≤ 20).

Errors (if any).

2. Variant Table: Outputs a table with details for each variant where QUAL > 20: Chromosome, Position, Ref/Alt bases, Variant Type (Non-coding, Synonymous, Non-synonymous), Transcript ID, Protein location, Ref/Alt amino acids.

3. Bar Plot: Displays the proportions of non-coding, synonymous, and non-synonymous variants for variants with QUAL > 20.
