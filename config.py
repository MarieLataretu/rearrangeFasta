OUTPUT_DIR = '/mnt/dessertlocal/kons/Aspergillus_fumigatus_A1163/rearr_spades_n157/'#'/mnt/prostlocal/marie/chlamydia_comparison/test/'#chlamydia_psittaci/'

BLAST_OUTPUT_DIR = ''.join([OUTPUT_DIR, 'blast/'])
POI_OUTPUT_DIR = ''.join([OUTPUT_DIR, 'POI/'])

BLAST_BIN_PATH = '/mnt/prostlocal/programs/blast/ncbi-blast-2.4.0+/bin/'

### Genabnk search options
## sequence feature type (CDS, gene, mRNA, ...)
seq_feature_type = 'CDS'
## sequence feature qualifier (product, protein_id, locus_tag, db_xref, ...)
seq_feature_qualifier = 'product'