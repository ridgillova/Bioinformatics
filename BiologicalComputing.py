# import necessary modules
import argparse
import logging
import gffutils
import vcf
import os
import sqlite3
#import pyfaidx
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
#from Bio.SeqRecord import SeqRecord
#from Bio import SeqIO
import math
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

# setting up logging
logger = logging.getLogger()
# setting logging level
logger.setLevel(logging.INFO)
# creating handler and setting format for handler, adding handler to logger
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(handler)

# setting up argparse
parser = argparse.ArgumentParser(description = '')
parser.add_argument('--vcf', required = True, help = 'The name of the .vcf file containing information about variants.')
parser.add_argument('--gff', required = True, help = 'The name of the .gff file containing information about the genes.')
parser.add_argument('--genome', required = True, help = 'The name of the .fasta file containing the sequences of the chromosomes.')
parser.add_argument('-outP','--outputPrefix', default = '2479057', help = 'The desired prefix for all output files. (Default: 2479057).')
args = parser.parse_args()

# creating a function to change a base to its complement
def complement(base):
    base = str(base).upper()
    base = base.replace('A','t')
    base = base.replace('T','a')
    base = base.replace('C','g')
    base = base.replace('G','c')
    base = base.upper()
    return base

# assertions for complement() function
assert complement('A') == 'T' # checking one of the bases 
assert complement('a') == 'T' # checking a lower-case letter will also result in the correct complement (in upper case)
# checking other bases
assert complement('T') == 'A'
assert complement('C') == 'G'
assert complement('G') == 'C'
# checking a sequence (more than one base) as well as the other bases
assert complement('ATCG') == ('TAGC')

# defining input files from args
# checking whether file can be found, and has expected suffix (specific file format), otherwise inform user and exit
if os.path.isfile(args.vcf):
    if args.vcf.endswith(('.vcf','.vcf.gz')):
        vcf_file = args.vcf # 'assessmentData.vcf.gz'
    else:
        logger.error(f'The input file {args.vcf} is not in vcf format. Please check and try again. (Accepted suffix: .vcf/.vcf.gz).')
        raise SystemExit(1)
else:
    logger.error(f'File {args.vcf} not found. Please check and try again.')
    raise SystemExit(1)
if os.path.isfile(args.gff):
    if args.gff.endswith('.gff'):
        gff_file = args.gff # 'PlasmoDB-54_Pfalciparum3D7.gff'
    else:
        logger.error(f'The input file {args.gff} is not in gff format. Please check and try again. (Accepted suffix: .gff).')
        raise SystemExit(1)
else:
    logger.error(f'File {args.gff} not found. Please check and try again.')
    raise SystemExit(1)
if os.path.isfile(args.genome):
    if args.genome.endswith(('.fasta', '.fa')):
        genome_file = args.genome # 'PlasmoDB-54_Pfalciparum3D7_Genome.fasta'
    else:
        logger.error(f'The input file {args.genome} is not in fasta format. Please check and try again. (Accepted suffix: .fasta).')
        raise SystemExit(1)
else:
    logger.error(f'File {args.genome} not found. Please check and try again.')
    raise SystemExit(1)
# add appropriate suffix to output files
output = str(args.outputPrefix) + '.tsv' # '2479057.tsv'
output_plot = str(args.outputPrefix) + '.png' #'2479057_plot.png'

# replacing .gff suffix by .db to create the database (or check whether it already exists)
db = gff_file.replace('gff', 'db')

# information for user about input files and output files
logger.info(f'The following input files are being used:\nvcf_file = {vcf_file}\ngff_file = {gff_file}\ngenome_file = {genome_file}')
logger.info(f'The following output files will be created:\nTsv file = {output}\nBar plot = {output_plot}')

# creating database from gff file
if not os.path.isfile(db): # check if database exists, if not create database
    logger.info(f'Creating database {db}.') # inform user that database is being created
    try:
        db = gffutils.create_db(gff_file, dbfn = db, force = True, keep_order = True) # delete force = True
    except sqlite3.OperationalError: # inform user if database cannot be created 
        logger.error(f'Cannot create database {db} from file {gff_file}.')
        raise SystemExit(1)
else: # if already exists, connect to database and inform user
    logger.info(f'Connecting to existing database {db}.')
    try:
        db = gffutils.FeatureDB(db, keep_order = True)
    except sqlite3.OperationalError:
        logger.error(f'Cannot connect to database {db}. Please check and try again.')
        raise SystemExit(1)

# read vcf file to be able to iterate through
try:
    vcf_file = vcf.Reader(filename = vcf_file)
except FileNotFoundError: # if file not found, inform user and exit
    logger.error(f'The input file {vcf_file} not found. Please check and try again.')
    raise SystemExit(1)

# creating counts to be able to inform user and make bar plot
count_variants = 0
count_qual_above = 0
count_qual_below = 0
count_coding = 0
count_noncoding = 0
count_synonymous = 0
count_nonsynonymous = 0
# open output file to be able to write in
with open(output, 'wt') as out:
    out.write('Chrom\tPos\tRef\tAlt\tType\tTranscript\tProtein Location\tRef AA\tAlt AA\n') # writing the header
    # for each record in vcf_file
    for record in vcf_file:
        count_variants += 1 # count variants processed to be able to inform user
        # check whether quality is above 20
        try:
            if record.QUAL > 20:
                count_qual_above += 1 # count all variants with QUAL > 20
                # check whether variant is in a coding region - if list contains more than one CDS
                if len(list(db.region(seqid = record.CHROM, start = record.POS, end = record.POS, featuretype = 'CDS'))) > 0:
                    count_coding += 1 # count variants in coding regions
                    # get the CDS where the variant is located 
                    for feature in db.region(seqid = record.CHROM, start = record.POS, end = record.POS, featuretype = 'CDS'):
                        # get the parent (transcript) of the CDS features
                        for parent in db.parents(feature.id, featuretype = 'mRNA'):
                            seq = '' # for each transcript create an empty string to store the sequence of the whole mature DNA containing the variant
                            snp_coordinate = 0 # create a variable to ge the SNP coordinate relative to start of the mature DNA
                            # if the transcript is on the positive strand:
                            if parent.strand == '+':
                                # get all the CDS regions corresponding to the given transcript (that create the mature DNA)
                                for child in db.children(parent.id, featuretype = 'CDS', order_by = 'start'):
                                    seq += child.sequence(genome_file, use_strand = True) # add the sequence corresponding to that CDS using the genome_file
                                    # if the CDS does not contain the variant and occurs before the variant 
                                    if child.start < record.POS and child.end < record.POS:
                                        # add the length of the CDS to the snp_coordinate
                                        snp_coordinate += (child.end - child.start + 1) # adding 1 to the difference between end and start of the CDS to get the length of the CDS
                                    # if the CDS contains the variant
                                    elif child.start < record.POS and child.end > record.POS:
                                        # add the number of bases up to the position of the variant to the snp_coordinate
                                        snp_coordinate += (record.POS - child.start)
                            # if the transcript is located on the negative strand, get the sequence and SNP coordinate in a similar manner as for the positive strand
                            # but using the reverse complement of the sequence from the genome_file
                            # and count the coordinate from the 5' end (i.e. from the child.end)
                            if parent.strand == '-':
                                for child in db.children(parent.id, featuretype = 'CDS', order_by = 'start', reverse = True):
                                    seq += child.sequence(genome_file, use_strand = True) # setting use_strand = True creates the reverse complement of the sequence as its on the negative strand
                                    if child.start > record.POS and child.end > record.POS: # adding to the SNP coordinate if the variant position is lower than the CDS range
                                        snp_coordinate += (child.end - child.start + 1)
                                    elif child.end > record.POS and child.start < record.POS: # if the CDS contains the SNP
                                        snp_coordinate += (child.end - record.POS) # add the number of bases from the end up to the SNP position
                            # create a MutableSeq object - to be able to change the reference base to the alternative base
                            mutableSeq = MutableSeq(seq)
                            # change the base at the snp_coordinate to the alternative base (use the complement function for the negative strand)
                            if parent.strand == '+':
                                mutableSeq[snp_coordinate] = str(record.ALT[0]) # have to index record.ALT as it's a list and convert to a string
                            if parent.strand == '-':
                                mutableSeq[snp_coordinate] = complement(str(record.ALT[0]))
                            # save as a Seq object to be able to translate the sequence (for both the reference and mutated sequence)
                            seq = Seq(seq).translate()
                            seq_mutated = Seq(mutableSeq).translate()
                            # check whether the protein is valid
                            try:
                                assert seq.startswith('M') and seq.endswith('*') and seq.count('*') == 1
                            # if not valid, inform the user and continue 
                            except AssertionError:
                                logger.warning(f'The protein sequence for {parent.id} is not valid.')
                                continue
                            # calculate the position of the amino acid containing the variant (substract 1 from the snp_coordinate so that multiples of 3 corresponds to the first position of a codon)
                            # the resulting aa_position will correspond to python-indexing (i.e., 0 corresponds to the first amino acid, 1 to the second etc.)
                            aa_position = math.ceil((snp_coordinate-1)/3) # corresponds to python-indexing
                            # get the amino acid at the mutation position (before and after mutation)
                            aa = seq[aa_position]
                            aa_mutated = seq_mutated[aa_position]
                            # if the amino acids are the same before and after mutation 
                            if aa == aa_mutated:
                                count_synonymous += 1 # count synonymous variants
                                type = 'synonymous' # set type to 'synonymous'
                                aa_mutated = 'NA' # set aa_mutated to NA 
                            # else (amino acids not the same), count non-synonymous variants and set type to 'non-synonymous'
                            else:
                                count_nonsynonymous += 1
                                type = 'non-synonymous'
                            # description = '{record.CHROM}: {feature.start} - {feature.end}'
                            # write all necessary info about variant into the output file
                            out.write(f'{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\t{type}\t{parent.id}\t{aa_position+1}\t{aa}\t{aa_mutated}\n')
                # if no CDS regions correspond to the variant
                else:
                    count_noncoding += 1 # count non-coding variants
                    type = 'non-coding' # set type as 'non-coding'
                    # write variant out into output file (setting information about transcript/protein as NA)
                    out.write(f'{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\t{type}\tNA\tNA\tNA\tNA\n')       
            # count variants that have QUAL <= 20
            else:
                count_qual_below += 1
        except TypeError:
            logger.error(f'File {vcf_file} has unexpected values in QUAL column. Please check file and try again.')

count_total = count_nonsynonymous + count_synonymous + count_noncoding
prop_nonsynonymous = count_nonsynonymous/count_total
prop_synonymous = count_synonymous/count_total
prop_noncoding = count_noncoding/count_total

# create a data frame contaning the counts for the different types of variants
type_counts = pd.DataFrame({'Type': ['synonymous', 'non-synonymous', 'non-coding'],
                            'Proportion': [prop_synonymous, prop_nonsynonymous, prop_noncoding]})

# create a bar plot, save figure and clear figure
sns.barplot(data = type_counts, x = 'Type', y = 'Proportion') # vertical bar plot
plt.savefig(output_plot) # using the name set by user or default from argparse
plt.clf()

# logging - info
logger.info(f'{count_variants} variants have been processed.')
logger.info(f'{count_qual_below} variants have QUAL <= 20.')
logger.info(f'{count_qual_above} variants have QUAL > 20.')
logger.info(f'There are {count_coding} variants located within coding regions, of which {count_synonymous} are synonymous and {count_nonsynonymous} are non-synonymous.')
logger.info(f'There are {count_noncoding} variants located within non-coding regions (including pseudogenes).')
logger.info(f'The output tab-separated file containing information about the variants is in {output}.')
logger.info(f'The output bar plot showing the proportion of synonymous, non-synonmous and non-coding variants is in {output_plot}.')
