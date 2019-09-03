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


def find_roi_in_gff_file(anno_with_POI_file, product_of_interest, output_file=None, fasta_file=None):
    breakpoint = None
    reverse_complement = False
    return((int(breakpoint), reverse_complement))

def find_roi_in_genbank_file(anno_with_POI_file, product_of_interest, output_file=None):
    target_id = None
    breakpoint = None
    reverse_complement = False

    ### search for the product of interest in the annotation file
    with open(anno_with_POI_file, 'r') as gbk:
        for seq_record in SeqIO.parse(gbk, 'genbank'):
            
            numb_of_ROIs = 0
            for seq_feature in seq_record.features:
                if seq_feature.type == config.seq_feature_type:
                    if config.seq_feature_qualifier in seq_feature.qualifiers:
                        if product_of_interest in seq_feature.qualifiers[config.seq_feature_qualifier]:
                            numb_of_ROIs += 1

                            if numb_of_ROIs == 1:
                                breakpoint = seq_feature.location.start
                                target_id = seq_record.id

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
    return(str(target_id), (int(breakpoint), reverse_complement))

def find_roi_in_fasta_file(fasta_file, roi_fasta, output_blast_name=None):
    reverse_complement = False

    ### read in the genome file
    num_of_records = 0
    with open(fasta_file) as fas:
        for record in SeqIO.parse(fas, 'fasta'):
            num_of_records += 1

    ### read in the roi file
    num_of_records = 0
    with open(roi_fasta) as fas:
        for record in SeqIO.parse(fas, 'fasta'):
            num_of_records += 1
            roi_id = record.id
    assert num_of_records == 1, f'The file {roi_fasta} is a multiple fasta file. Please use a single fasta file with only one sequence.'

    if not (os.path.isfile(fasta_file + '.nhr') and os.path.isfile(fasta_file + '.nsq') and os.path.isfile(fasta_file + '.nin')):
        print ('start makeblastdb')
        makeblastdbStr = f'makeblastdb -dbtype nucl -in {fasta_file}'
        os.system(makeblastdbStr)
        print ('end makeblastdb')

    if output_blast_name:
        # make dir structure and check if file exist
        os.makedirs(os.path.dirname(output_blast_name), exist_ok=True)
        if os.path.isfile(output_blast_name):
            print(f'WARNING: overriding {output_blast_name}')
    else:
        output_blast_name = os.path.join(os.path.dirname(os.path.abspath(fasta_file)), f'{os.path.splitext(fasta_file)[0]}_vs_{roi_id}_blast.tsv')


    print ('start blast')
    blastStr = f'blastn -out {output_blast_name} -outfmt \'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\' -query {roi_fasta} -db {fasta_file} -evalue 1e-10'
    os.system(blastStr)
    print ('end blast')

    ## read in the blast result
    with open(output_blast_name) as csvf:
        table = [row for row in csv.reader(csvf, delimiter='\t')]
    if len(table) == 1:

        ## check the stand of the product of interest
        if int(table[0][8]) > int(table[0][9]):
            return((str(table[0][0]), int(table[0][8])-1, True))
        else:
            return((str(table[0][0]), int(table[0][8])-1, False))

    elif len(table) == 0:
        print ('no blast hit')
        sys.exit(0)
    else:
        print ('end find product of interest: found more than one blast hit; not implemented yet')
        sys.exit(0)