### Genabnk search options
## sequence feature type (CDS, gene, mRNA, ...)
seq_feature_type = 'CDS'
## sequence feature qualifier (product, protein_id, locus_tag, db_xref, ...)
seq_feature_qualifier = 'product'


### blast parameters
# formtat is set to "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"; needed for further processing
makeblastdb_args = '-dbtype nucl'
blastn_args = '-evalue 1e-10'