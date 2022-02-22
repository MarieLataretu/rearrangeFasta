import sys
import os
import csv
import config
import warnings
import subprocess
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
    numb_of_ROIs = 0

    ### search for the product of interest in the annotation file
    with open(anno_with_POI_file, 'r') as gbk:
        for seq_record in SeqIO.parse(gbk, 'genbank'):
            
            for seq_feature in seq_record.features:
                if seq_feature.type == config.seq_feature_type:
                    if config.seq_feature_qualifier in seq_feature.qualifiers:
                        if product_of_interest in list(map(lambda x: x.lower(), seq_feature.qualifiers[config.seq_feature_qualifier])):
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
    assert breakpoint != None, f"Ooops, something went wrong. No breakpoint was found."
    return(str(target_id), (int(breakpoint), reverse_complement))

def find_roi_in_fasta_file(fasta_file, roi_fasta, output_blast_name=None):
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
    
    makeblastdb_log_file = f'{fasta_file}.makeblastdb.log'
    if not (os.path.isfile(fasta_file + '.nhr') and os.path.isfile(fasta_file + '.nsq') and os.path.isfile(fasta_file + '.nin')):
        print ('start makeblastdb')
        makeblastdb_cmd = f'makeblastdb -in {fasta_file} {config.makeblastdb_args}'
        with open(makeblastdb_log_file, 'w') as log:
            returncode_makeblastdb = subprocess.run(makeblastdb_cmd, shell=True, stdout=log, stderr=log).returncode
        assert returncode_makeblastdb == 0, f'ERROR: makeblastdb failed with exit code {returncode_makeblastdb}.'
        print ('end makeblastdb')

    if output_blast_name:
        # make dir structure and check if file exist
        os.makedirs(os.path.dirname(output_blast_name), exist_ok=True)
        if os.path.isfile(output_blast_name):
            print(f'WARNING: overriding {output_blast_name}')
    else:
        output_blast_name = os.path.join(os.path.dirname(os.path.abspath(fasta_file)), f'{os.path.splitext(fasta_file)[0]}_vs_{roi_id}_blast.tsv')
    blastn_log_file = f'{output_blast_name}.blastn.log'

    print ('start blast')
    blast_cmd = f'blastn -out {output_blast_name} -query {roi_fasta} -db {fasta_file} {config.blastn_args} -outfmt \'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\''
    blastn = subprocess.run(blast_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if blastn.stdout or blastn.stderr:
        with open(blastn_log_file, 'wb') as log:
            log.write(blastn.stdout)
            log.write(blastn.stderr)
    assert blastn.returncode == 0, f'ERROR: blastn failed with exit code {blastn.returncode}.'
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