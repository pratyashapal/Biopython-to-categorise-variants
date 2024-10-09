# Importing necessary modules
import vcf
import gffutils
import os
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
import logging
import argparse
import sqlite3

# Parsing the arguments
parser = argparse.ArgumentParser(description='Process VCF data.')
parser.add_argument('--vcf', dest='vcf_filename', required=True, help='Input VCF filename')
parser.add_argument('--gff', dest='gff_filename', required=True, help='Input GFF filename')
parser.add_argument('--fasta', dest='fasta_filename', required=True, help='Input fasta filename')
args = parser.parse_args()

# Getting the root logger instance
# Using loggerinfo to print the necessary information about the data 
# Setting the logging level to INFO
# Creating a FileHandler to handle the logging output to a file 
# Creating a Formatter to specify the format of the log messages
# Adding the FileHandler with the specified Formatter to the logger
loggerinfo = logging.getLogger('loggerinfo')
loggerinfo.setLevel(logging.INFO) 
ch = logging.FileHandler('2875659_log_file.log') 
ch.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s')) 
loggerinfo.addHandler(ch)

# Using loggererror to print the error messages
loggerError = logging.getLogger('loggerError')
loggerError.setLevel(logging.ERROR)
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
loggerError.addHandler(sh)

# Logging the VCF. GFF, Fasta files used
loggerinfo.info(f'The VCF file provided is:  {args.vcf_filename}\n')
loggerinfo.info(f'The GFF file provided is:  {args.gff_filename}\n')
loggerinfo.info(f'The Fasta file provided is: {args.fasta_filename}\n')

# Reading the VCF file
try: 
    vcfReader = vcf.Reader(filename=args.vcf_filename)
# Handling the error if the file is wrong and raising system exit
except FileNotFoundError:
    loggerError.error(f'Error reading VCF file: {args.vcf_filename}. Please check and try again\n')
    raise SystemExit(1)

# Creating the database file name by removing '.gff' and adding '.db' to the gff file
db_filename = os.path.splitext(args.gff_filename)[0] + ".db"

# Creating a new database if it doesn't exist
if not os.path.isfile(db_filename):
    try:
    #print('Creating database...\n')
        db = gffutils.create_db(args.gff_filename, dbfn=db_filename, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    except sqlite3.OperationalError:
        # Handling and logging errors
        loggerError.error(f'Cannot create database {db_filename}\n')
        raise SystemExit(1)
    except ValueError:
        loggerError.error(f'Cannot create database {db_filename}\n')
        raise SystemExit(1)
    
# If database already present, connecting to it
else:
    #print('Connecting to existing database.\n')
    try:
        db = gffutils.FeatureDB(db_filename, keep_order=True)
    except ValueError:
        loggerError.error(f'Database {db_filename} could not be read\n')
        raise SystemExit(1)
        

# Creating functions: 
# Getting the location till the cds with variant is met 
def get_sequence_and_location(child,args,seq, loc, variant):
    """
    Retrieve the sequence and location information from the child CDS feature.

    Parameters:
    - child: The current CDS feature.
    - variant: The variant feature.
    - args: Command-line arguments containing the fasta filename.

    Returns:
    - seq: Concatenated CDS sequence.
    - loc: Concatenated location information.
    """

    seq += child.sequence(args, use_strand=True)

    if child == variant:
        return seq, loc, True
    else:

        loc += child.sequence(args, use_strand=True)

        return seq, loc, False

def get_locations(mutated_position, translated_seq, new_translated_seq):
    """
    Retrieve the reference and alternative amino acids based on the protein location.

    Parameters:
    - mutated_position : position of the variant
    - translated_seq: Original translated sequence.
    - new_translated_seq: Mutated translated sequence.

    Returns:
    - ref_aa: Reference amino acid.
    - alt_aa: Alternative amino acid.
    - protein_location : location of variant in translated sequence
    """
    ref_aa = 'NA'
    alt_aa = 'NA'
    '''
    The mutated_position likely denotes a location within 
    the coding sequence influenced by the mutation. 
    Since amino acids are commonly indexed 
    starting from 1. The division by 3 reflects the idea that
    each codon, composed of three nucleotides, corresponds to
    a single amino acid. The addition of 1 is employed to facilitate 
    the transition from 0-based indexing to 
    1-based indexing.
    '''
    protein_loc = int((mutated_position) // 3 + 1)

    if 0 <= protein_loc <= len(translated_seq) and 0 <= protein_loc <= len(new_translated_seq):
        # Reverting the 1 base indexed protein location to 0 base indexing by substracting 1
        # Obtaining the Reference and Alternate alleles
        ref_aa = translated_seq[protein_loc - 1]
        #print(ref_aa)
        alt_aa = new_translated_seq[protein_loc - 1]
        #print(alt_aa)

    return ref_aa, alt_aa, protein_loc

# Function to adjust the sequence length to be a multiple of three
def adjust_sequence_length(sequence):
    # Calculating the number of extra bases
    extra_bases = len(sequence) % 3
    # Trimming the sequence if it's not a multiple of three
    if extra_bases != 0:
        sequence = sequence[:-extra_bases]
    return sequence

# Function to translate a DNA sequence to a protein sequence
def translate_sequence(seq):
    seq = adjust_sequence_length(seq)
    # Translating the sequence to a protein sequence
    translated_seq = Seq(seq).translate()
    return translated_seq

# Naming the output files
output_file = '2875659_BCPY_OUTPUT'
bar_plot_filename = '2875659_bar_plot.png'

# Creating a dictionary to store the number of non coding, synonymous and non synonymous variants
mutation_counts = {'non coding': 0, 'synonymous': 0, 'non synonymous': 0}

# Initialising the variables to be used further
low_qual = 0
trans = 'NA'
protein_loc = 'NA'
translated_seq = 'NA'
new_translated_seq = 'NA'
ref_aa = 'NA'
alt_aa = 'NA'

try:
    # Opening the output file in write mode ('w'), using a context manager
    with open(output_file, 'w') as f:
        
        # Writing a header to the output file
        f.write("CHROM\tPOS\tREF\tALT\tType\tTranscript\tProtein Location\tRef AA\tAlt AA\n")
        
        # Iterate over each variant record in the VCF file
        for record in vcfReader:
            
            # Checking if the quality of the variant is less than or equal to 20 (QUAL field in the VCF). 
            # This will be recorded in the log file
            if record.QUAL <= 20:
                low_qual += 1
                
            # Processing variants with quality greater than 20
            if record.QUAL > 20:
                mutation_type = 'non coding'
                
                # Adding the non coding counts to the mutation count dictionary
                mutation_counts['non coding'] += 1
                
                # Retrieving coding sequences overlapping with the variant
                coding = list(db.region(seqid=record.CHROM, start=record.POS, end=record.POS, featuretype='CDS'))
                
                # Format for the coding sequences
                if coding:
                    for variant in coding:
                        
                        # Getting the cds parents (going up the hierarchy)
                        for transcript in db.parents(variant, featuretype='mRNA'):
                            seq = ''
                            loc = ""
                            trans = transcript.id
                            
                            # Determining the strand of the transcript
                            if transcript.strand == '+':
                                
                                # For positive strand, concatenate CDS sequences in order
                                for child in db.children(transcript, featuretype='CDS', order_by='start'):
                                    
                                    # Calling the function to get the appropriate concatenated sequence.
                                    # Obtaining the relative location of the variant (loc)
                                    # cds_true here implies True in the function i.e if the cds is the cds with the variant
                                    seq, loc, cds_true =get_sequence_and_location(child,args.fasta_filename,seq,loc,variant)
                                    if cds_true:
                                        break
                                        
                                if seq:
                                    
                                    start_pos = variant.start
                                    
                                    # Calculating the position of the variant
                                    mutated_position = record.POS - start_pos + len(loc)
                                    
                                    # Creating a mutable sequence for modification
                                    mutableseq = MutableSeq(seq)
                                    
                                    # Iterating through all the alternate alleles in record
                                    for alt in record.ALT:
                                        mutableseq[mutated_position] = str(alt)
                                        newseq = str(mutableseq)
                                        
                                    # Translating the original and mutated sequences
                                    translated_seq = translate_sequence(seq)
                                    new_translated_seq = translate_sequence(newseq)
                                    
                                    # Extracting the reference, alterante allele and the protein location
                                    ref_aa, alt_aa,protein_loc = get_locations(mutated_position, translated_seq, new_translated_seq)

                            elif transcript.strand == '-':
                                
                                # For negative strand, concatenating CDS sequences in reverse order
                                trans = transcript.id
                                for child in db.children(transcript, featuretype='CDS', order_by='start', reverse=True):
                                    seq, loc, cds =get_sequence_and_location(child,args.fasta_filename,seq,loc,variant)
                                    if cds:
                                        break
                                        
                                if seq:
                                    end_pos = variant.end                                    
                                    mutated_position =  (end_pos - record.POS + len(loc)) # Calculating the position of the variant
                                    mutableseq = MutableSeq(seq) # Creating a mutable sequence for modification
                                    
                                    # Iterating through all the alternate alleles in record
                                    for alt in record.ALT:
                                        
                                        # Getting the complement of the allele since they are in respect to the + strand
                                        substituted_seq = Seq(alt.sequence).complement()
                                        mutableseq[mutated_position ] = str(substituted_seq)
                                        newseq = str(mutableseq)
                                        
                                        # Translating the original and mutated sequences
                                        translated_seq = translate_sequence(seq)
                                        new_translated_seq = translate_sequence(newseq)
                                        
                                        # Extracting the reference, alterante allele and the protein location
                                        ref_aa, alt_aa,protein_loc = get_locations(mutated_position, translated_seq, new_translated_seq)
                            
                            # Handling the case where the the transcript is neither '+' nor '-' 
                            # Malformatted file
                            else:
                                # Logging the error
                                loggerError.error(f'GFF file {args.gff_filename} is malformatted')
                                raise SystemExit()
                            
                            # Determining mutation type based on translated sequences
                            if new_translated_seq == translated_seq:
                                mutation_type = 'synonymous'
                                # Adding the synonymous counts to the mutation count dictionary
                                mutation_counts['synonymous'] += 1
                            else:
                                mutation_type = 'non synonymous'
                                # Adding the non synonymous counts to the mutation count dictionary
                                mutation_counts['non synonymous'] += 1
                                
                        # Writing a tab separated file with the desired values
                        f.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\t{mutation_type}\t{trans}\t{protein_loc}\t{ref_aa}\t{alt_aa}\n")
                        
                # Format for the non coding sequences
                else:
                    f.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\t{mutation_type}\tNA\tNA\tNA\tNA\n")
                    
    loggerinfo.info(f'Successful execution of the script. Writing to output file: {output_file}\n')
    
# Handling the error if the file cannot be opened and exiting
except FileNotFoundError:
    loggerError.error(f' Output File {output_file} cannot be opened for reading. Please check and try again\n')
    raise SystemExit(1)

# Logging the variants with count <=20
loggerinfo.info(f'The count of variants with QUAL <= 20: {low_qual}\n')

# Extracting mutation types and corresponding counts for plotting
# List of mutation types (e.g., 'non coding', 'synonymous', 'non synonymous')
mutation_types = list(mutation_counts.keys()) 

# List of counts corresponding to each mutation type
counts = list(mutation_counts.values()) 

# Creating a bar plot using Matplotlib
plt.bar(mutation_types, counts, color=['lightcoral', 'lime', 'lightseagreen']) # Customising the plot appearance
plt.xlabel('Mutation Type') # Setting label for x axis
plt.ylabel('Count') # Setting label for y axis
plt.title('Proportion of Variants with QUAL > 20') # Setting the title of the plot
plt.show() # Displaying the plot
plt.savefig(bar_plot_filename) # Saving the plot as a PNG image file

# Log information about writing to the output file and saving the bar plot
loggerinfo.info(f'Bar plot saved as: {bar_plot_filename}\n')
