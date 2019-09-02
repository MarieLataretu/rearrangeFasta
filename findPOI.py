import sys
import os
import csv
import config
import warnings
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature as sf

class NoMatchFoundError(Exception):
    def __init__(self, inp_file, target):
        self.inp_file = inp_file
        self.target = target
    def __str__(self):
        return (f'ERROR: No match of "{self.target}" found in file {self.inp_file}')

class NoMatchFoundInGenbankError(NoMatchFoundError):
    def __init__(self, inp_file, target, feature, quali):
        self.inp_file = inp_file
        self.target = target
        self.feature = feature
        self.quali = quali
    def __str__(self):
        return (f'ERROR: No match of "{self.target}" found in file {self.inp_file} with "{self.feature}" as feature type and "{self.quali}" as qualifier.\n To change the feature type or qualifier, please change the config.')

class MultipleMatchesFoundWarning(UserWarning):
    def __init__(self, inp_file, target):
        self.inp_file = inp_file
        self.target = target
    def __str__(self):
        return (f'Multiple matches of "{self.target}" found in file {self.inp_file}')


def find_poi_in_gff_file(anno_with_POI_file, product_of_interest, output_file=None, fasta_file=None):
    breakpoint = None
    reverse_complement = False
    return((int(breakpoint), reverse_complement))

def find_poi_in_genbank_file(anno_with_POI_file, product_of_interest, output_file=None):
    breakpoint = None
    reverse_complement = False

    ### search for the product of interest in the annotation file
    with open(anno_with_POI_file, 'r') as gbk:
        for seq_record in SeqIO.parse(gbk, 'genbank'):
            
            numb_of_ROIs= 0
            for seq_feature in seq_record.features:
                if seq_feature.type == config.seq_feature_type:
                    if config.seq_feature_qualifier in seq_feature.qualifiers:
                        if product_of_interest in seq_feature.qualifiers[config.seq_feature_qualifier]:
                            numb_of_ROIs += 1

                            if numb_of_ROIs == 1:
                                breakpoint = seq_feature.location.start

                                ## check the strand
                                if seq_feature.location.strand == -1:
                                    reverse_complement = True

                                if output_file:
                                    record_description = seq_record.description.replace(' ', '_')
                                    ## save the nucleotide sequence (e.g. to blast it against an another genome)
                                    if type(seq_feature.location) is sf.FeatureLocation:
                                        ## file to save the nucleotide sequence of the product of interest
                                        if reverse_complement:
                                            new_description = f"{record_description}-{'revcomp'}-:{seq_feature.location.start + 1}-{seq_feature.location.end}"

                                            new_nc_record = SeqRecord(
                                                seq_record.seq[seq_feature.location.start:seq_feature.location.end].reverse_complement(),
                                                id=product_of_interest, description=new_description)
                                        else:
                                            new_description = f"{record_description}:{seq_feature.location.start + 1}-{seq_feature.location.end}"
                                            new_nc_record = SeqRecord(seq_record.seq[seq_feature.location.start:seq_feature.location.end], id=product_of_interest, description=new_description)

                                        # make dir structure and check if file exist
                                        os.makedirs(os.path.dirname(output_file), exist_ok=True)
                                        if os.path.isfile(output_file):
                                            print(f'WARNING: overriding {output_file}')
                                        SeqIO.write(new_nc_record, output_file, 'fasta')
                                    else:
                                        print('product of interest has a CompoundLocation; this is not implemented yet')
                                        sys.exit(1)
    if numb_of_ROIs == 0:
        raise NoMatchFoundInGenbankError(anno_with_POI_file, product_of_interest, config.seq_feature_type, config.seq_feature_qualifier)
    elif numb_of_ROIs > 1:
        warnings.warn(MultipleMatchesFoundWarning(anno_with_POI_file, product_of_interest))
    assert breakpoint, f"Ooops, something went wrong. No breakpoint was found."
    return((int(breakpoint), reverse_complement))

def find_poi_in_fasta_file(genome_fasta_file, seqeunce_of_interest_fasta_file, output_blast_name):
    reverse_complement = False
    save_contig = False

    if (os.path.isdir(config.BLAST_OUTPUT_DIR)):
        os.mkdir(config.BLAST_OUTPUT_DIR)
        print('blast output directory created')

    ### read in the genome file
    num_of_records = 0
    with open(genome_fasta_file) as fas:
        for record in SeqIO.parse(fas, 'fasta'):
            if ('plasmid' not in record.description):
                genome = record
                num_of_records = num_of_records + 1
    if (num_of_records > 1):
        save_contig = True

    if (not os.path.exists(config.BLAST_OUTPUT_DIR)):
        os.makedirs(config.BLAST_OUTPUT_DIR)
    # out_blast_file = config.BLAST_OUTPUT_DIR + strain_name_with_POI + '_poi_vs_' + strain_name_to_search_in + '_genome.tab'
    out_blast_file = config.BLAST_OUTPUT_DIR + output_blast_name +'.tab'

    if (not (os.path.isfile(genome_fasta_file + '.nhr') and os.path.isfile(genome_fasta_file + '.nsq') and os.path.isfile(genome_fasta_file + '.nin'))):
        print ('start makeblastdb')
        makeblastdbStr = config.BLAST_BIN_PATH + 'makeblastdb -dbtype nucl -in ' + genome_fasta_file
        os.system(makeblastdbStr)
        print ('end makeblastdb')

    print ('start blast')
    blastStr = config.BLAST_BIN_PATH + 'blastn -out ' + out_blast_file + ' -outfmt \'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\' -query ' + seqeunce_of_interest_fasta_file + ' -db ' + genome_fasta_file + ' -evalue 1e-10'
    os.system(blastStr)
    print ('end blast')

    ## read in the blast result
    with open(out_blast_file) as csvf:
        table = [row for row in csv.reader(csvf, delimiter='\t')]
    if (len(table) is 1):

        ## check the stand of the product of interest
        if (int(table[0][8]) > int(table[0][9])):
            reverse_complement = True
        # return the breakpoint, if is is reverse complemented and on which contig the breakpoint is, if necessery
        if(save_contig):
            return ((int(table[0][8])-1, reverse_complement, table[0][1]))
        else:
            return ((int(table[0][8])-1, reverse_complement))
    elif(len(table) is 0):
        print ('no blast hit')
        sys.exit(0)
    else:
        print ('end find product of interest: found more than one blast hit; not implemented yet')
        sys.exit(0)