import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

genome_file = '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_gallinacea/14DC101_Cgallinacea_genome_rear.fasta'

genomes = []
with open(genome_file) as fas:
    for record in SeqIO.parse(fas, 'fasta'):
        genomes.append(record)

new_genome = SeqRecord(genomes[0].seq+genomes[2].seq+genomes[1].seq, id='14DC101_Cgallinacea', description='linked genome sequence')

SeqIO.write(new_genome, os.path.dirname(os.path.abspath(genome_file)) + '/14DC10_Cgallinacea_genome_rear_linked_132.fasta', 'fasta')