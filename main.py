import os
import sys
import argparse
from findPOI import NoMatchFoundError
from findPOI import MultipleMatchesFoundWarning
from findPOI import find_poi_in_gff_file
from findPOI import find_poi_in_genbank_file
from findPOI import find_poi_in_fasta_file
from rearranging import rearrange_fasta_file
from rearranging import rearrange_genbank_file
from rearranging import rearrange_gff_file
from rearranging import get_sequence_length_fasta
from rearranging import get_sequence_length_genbank

"""This script rearranges a genome (fasta file) and, if given, a annotaion (Genbank file) with respect to a product
of interest (POI).

Output paths can be specified in config.py. 

Usage:
python main.py "strain_name_with_POI" "genome_with_POI_file" "anno_with_POI_file" "product_of_interest"
["strain_name_to_search_in"] ["genome_to_search_in"] ["anno_to_search_in"]

Examples:
a) given a genome and annotation file of a species called 'fuubarus' and a annotated product 'muh' 
-> python main.py "fuubarus" "fuubarus.fasta" "fuubarus.gbf" "muh"

b) given a genome file of a species called 'barbarus' and a product 'muh', that is not annotated. take the genome and
annotation of 'fuubarus' to get the sequence of 'muh' and run blast to find it in 'barbarus'
-> python main.py "fuubarus" "fuubarus.fasta" "fuubarus.gbf" "muh" "barbarus" "barbarus.fasta"

c) given a genome file of a species called 'barbarus' and a product 'muh', that is not annotated namely in the 
annotation. 
take the genome and annotation of 'fuubarus' to get the sequence of 'muh' and run blast to find it in 'barbarus'.
-> python main.py "fuubarus" "fuubarus.fasta" "fuubarus.gbf" "muh" "barbarus" "barbarus.fasta" "barbarus.gbf"

Used tools: ncbi-blast-2.4.0+


Author: Marie Lataretu
E-Mail: marie.lataretu@uni-jena.de
"""
GENBANK_EXTENSIONS= ['.gb', '.gbff', '.genebank', '.gbk']
GFF_EXTENSIONS = ['.gff', '.gff3', '.gtf'] # works gtf?
FASTA_EXTENSIONS = ['.fasta', '.fa', '.fn']

class ExtensionNotFoundError(Exception):
    def __init__(self, in_file, gb_ext=GENBANK_EXTENSIONS, gff_ext=GFF_EXTENSIONS, fa_ext=FASTA_EXTENSIONS):
        self.in_file = in_file
        self.gb_ext = gb_ext
        self.gff_ext = gff_ext
        self.fa_ext = fa_ext
    def __str__(self):
        return (f"The file extension of {self.in_file} is unknown.\nAllowed extensions: {', '.join(self.gb_ext)} for GenBank files, {', '.join(self.gff_ext)} for GFF files and {', '.join(self.fa_ext)} Fasta files.")

def name_search(args):
    extension = os.path.splitext(args.annotation_file)[1].lower()

    searchResult = None
    try:
        if extension in GENBANK_EXTENSIONS:
            try:
                searchResult = find_poi_in_genbank_file(args.annotation_file, args.region_of_interest.replace('\'', ''), args.save_roi_sequence)
            except NoMatchFoundError as err:
                print(err)
                sys.exit(1)
        elif extension in GFF_EXTENSIONS:
            searchResult = find_poi_in_gff_file(args.annotation_file, args.region_of_interest.replace('\'', ''), args.save_roi_sequence)
        else:
            raise ExtensionNotFoundError(args.annotation_file)
    except ExtensionNotFoundError as err:
        print(err)
        sys.exit(1)

    print(searchResult)

def structural_search(args):
    searchResult = find_poi_in_fasta_file(args.genomeFile, args.sequenceOfInterest, args.blastOutName)

desc = 'This script rearranges a Fasta, Genbank or GFF file with respect to a anotated region of interest or specified sequence in Fasta format.'
# usage = 'python main.py "strain_name_with_POI" "genome_with_POI_file" "anno_with_POI_file" "product_of_interest" ["strain_name_to_search_in"] ["genome_to_search_in"] ["anno_to_search_in"]'
# usage = ' python main.py [-nameSearch|-structuralSearch]\nfor name search: python main.py -nameSearch -annotationFile file -POIName gene\nfor structural search: python main.py -structuralSearch'
# epilog = 'Output paths can be specified in config.py.\n\nExamples:\na) given a genome and annotation file of a species called 'fuubarus' and a annotated product 'muh'\n> python main.py "fuubarus" "fuubarus.fasta" "fuubarus.gbf" "muh"\n\nb) given a genome file of a species called 'barbarus' and a product 'muh', that is not annotated. take the genome and annotation of 'fuubarus' to get the sequence of 'muh' and run blast to find it in 'barbarus'\n-> python main.py "fuubarus" "fuubarus.fasta" "fuubarus.gbf" "muh" "barbarus" "barbarus.fasta"\n\nc) given a genome file of a species called 'barbarus' and a product 'muh', that is not annotated namely in the annotation. take the genome and annotation of 'fuubarus' to get the sequence of 'muh' and run blast to find it in 'barbarus'.\n-> python main.py "fuubarus" "fuubarus.fasta" "fuubarus.gbf" "muh" "barbarus" "barbarus.fasta" "barbarus.gbf"'
epilog = 'Author: Marie Lataretu\nE-Mail: marie.lataretu@uni-jena.de\nGitHub: https://github.com/MarieLataretu/rearrangeFasta'

parser = argparse.ArgumentParser(description=desc, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)

search_type_subparser = parser.add_subparsers(help='Search modi:')

name_search_parser = search_type_subparser.add_parser('nameSearch', help='search by an annotated feature')
name_search_parser.set_defaults(func=name_search)
required_name_search_args = name_search_parser.add_argument_group('Required arguments')
required_name_search_args.add_argument('--annotation_file', type=str, metavar='FILE', required=True, help='genbank or gff file')
required_name_search_args.add_argument('-roi', '--region_of_interest', type=str, metavar='STR', required=True, help='name of the annotated region of interest')
rearrange_target_name_search_args = name_search_parser.add_argument_group('Rearrange target file(s)', 'at least one file is required')
rearrange_target_name_search_args.add_argument('files', type=str, metavar='FILE', nargs='+', help='list of files to rearrange (possible file types: fasta, genbank, gff)')
output_name_search_options = name_search_parser.add_argument_group('Output options')
output_name_search_options.add_argument('--output_name', type=str, metavar='STR', help='name of the output file')
output_name_search_options.add_argument('--save_roi_sequence', type=str, metavar='FILE', help='name for the sequence outputfile of the annotated region of interest')

structural_search_parser = search_type_subparser.add_parser('structuralSearch', help='search by blasting a sequence')
structural_search_parser.set_defaults(func=structural_search)
structural_search_parser.add_argument('-annotationFile', type=str, metavar='FILE', help='genbank or gff file')
structural_search_parser.add_argument('-POIName', type=str, metavar='STR', help='name of the annotated region of interest')
structural_search_parser.add_argument('-sequenceOfInterest', type=str, metavar='FILE', help='nucleotide sequence for structural search')
structural_search_parser.add_argument('-blastOutName', type=str, metavar='STR', help='name for the outputfile of the structural search')
structural_search_parser.add_argument('-genomeFile', type=str, metavar='FILE', help='fasta file of the genome')
structural_search_parser.add_argument('-gffFile', type=str, metavar='FILE', help='gff file of the genome')
output_structural_search_options = structural_search_parser.add_argument_group('Output options')
output_structural_search_options.add_argument('-outputName', type=str, metavar='STR', help='name of the output file')
output_structural_search_options.add_argument('-savePOISequence', type=str, metavar='FILE', help='name for the sequence outputfile of the annotated region of interest')

rearrange_target = parser.add_argument_group('Target file to rearrange')
rearrange_target.add_argument('-rearrangeFasta', action='store_true', help='rearrange a fasta file')
rearrange_target.add_argument('-rearrangeGenbank', action='store_true', help='rearrange a genbank file')
rearrange_target.add_argument('-rearrangeGff', action='store_true', help='rearrange a gff file')

# args = parser.parse_args(['-nameSearch', '-annotationFile', 'TEST_Cavium_annotation_small.gbk', '-POIName', '\'delta-aminolevulinic acid dehydratase\''])
# args = parser.parse_args(['-nameSearch', '-annotationFile', '/mnt/prostlocal/marie/chlamydia_comparison/test/Cmu_Nigg.gbff',
#                           '-POIName', '\'delta-aminolevulinic acid dehydratase\'', '-savePOISequence', '-POIOutName', 'cmu_nigg_test'])
# args = parser.parse_args(['-structuralSearch', '-genomeFile', '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_genome.fasta',
#                           '-sequenceOfInterest', '/mnt/prostlocal/marie/chlamydia_comparison/test/POI/testseqpsittaci.fasta', '-blastOutName', 'blasttestpsittaci'])
# args = parser.parse_args(['-nameSearch', '-annotationFile', '/mnt/prostlocal/marie/chlamydia_comparison/test/Cmu_Nigg.gbff',
#                           '-POIName', '\'delta-aminolevulinic acid dehydratase\'',
#                             '-rearrangeFasta', '-genomeFile', '/mnt/prostlocal/marie/chlamydia_comparison/test/Cmu_Nigg.fasta',
#                           '-outputName', 'cmu_nigg_test'])
# args = parser.parse_args(['-structuralSearch', '-sequenceOfInterest', '/mnt/prostlocal/marie/chlamydia_comparison/test/POI/cmu_nigg_test.fasta',
#                           '-blastOutName', 'cmu_nigg_test',
#                             '-rearrangeFasta', '-genomeFile', '/mnt/prostlocal/marie/chlamydia_comparison/test/Cmu_Nigg.fasta',
#                           '-outputName', 'cmu_nigg_test'])
# args = parser.parse_args(['-nameSearch', '-annotationFile', '/mnt/prostlocal/marie/chlamydia_comparison/test/Cmu_Nigg.gbff',
#                           '-POIName', '\'delta-aminolevulinic acid dehydratase\'',
#                           '-rearrangeGenbank', '-outputName', 'testgenbank_name_psittaci'])
# args = parser.parse_args(['-nameSearch', '-annotationFile', '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_annotation.gbk',
#                           '-POIName', '\'delta-aminolevulinic acid dehydratase\'',
#                           '-rearrangeGff', '-gffFile', '/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_annotation.gff', '-outputName', 'testgff_name_psittaci'])
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
elif len(sys.argv) == 2:
    if sys.argv[-1] == 'nameSearch':
        name_search_parser.print_help(sys.stderr)
        sys.exit(1)
    elif sys.argv[-1] == 'structuralSearch':
        name_search_parser.print_help(sys.stderr)
        sys.exit(1)

args = parser.parse_args()
args.func(args)

# if(args.rearrangeFasta):
#     if(not(args.genomeFile or args.outputName)):
#         sys.exit('You need to define a fasta file with -genomeFile and a name for the output with -outputName.')
#     if(searchResult):
#         if(len(searchResult) is 3):
#             rearrange_fasta_file(args.genomeFile, searchResult[0], searchResult[1], args.outputName, searchResult[2])
#         else:
#             rearrange_fasta_file(args.genomeFile, searchResult[0], searchResult[1], args.outputName)
#     else:
#         sys.exit('Something went wrong with the search.')

# if(args.rearrangeGenbank):
#     if(not(args.annotationFile or args.outputName)):
#         sys.exit('You need to define a genbank file with -annotationFile and a name for the output with -outputName.')

#     if(searchResult):
#         rearrange_genbank_file(args.annotationFile, searchResult[0], searchResult[1], args.outputName)
#     else:
#         sys.exit('Something went wrong with the search.')

# if(args.rearrangeGff):
#     if(not (args.gffFile or args.outputName)):
#         sys.exit('You need to define a gff file with -gffFile and a name for the output with -outputName.')

#     if (searchResult):
#         if(args.genomeFile):
#             rearrange_gff_file(args.gffFile, get_sequence_length_fasta(args.genomeFile), searchResult[0], args.outputName)
#         else:
#             rearrange_gff_file(args.gffFile, get_sequence_length_genbank(args.annotationFile), searchResult[0], args.outputName)
#     else:
#         sys.exit('Something went wrong with the search.')

# def main(argv):
#     ### check and pares input parameters
#     if(len(argv) == 4 or len(argv) == 6 or len(argv) == 7):
#         strain_name_with_POI = argv[0].replace(' ', '_')
#         genome_with_POI_file = argv[1]
#         anno_with_POI_file = argv[2]
#         product_of_interest = argv[3]
#
#     else:
#         sys.exit('4, 6 or 7 arguments are expected')
#
#
#     if(len(argv) == 4):
#         ### search by name
#         breakpoint_and_reverse_complement = findPOI(strain_name_with_POI, genome_with_POI_file, anno_with_POI_file, product_of_interest)
#         rearrange_fasta_file(breakpoint_and_reverse_complement[0], breakpoint_and_reverse_complement[1], strain_name_with_POI, genome_with_POI_file, product_of_interest, anno_with_POI_file)
#
#     elif(len(argv) > 4):
#         ### structural search
#         strain_name_to_search_in = argv[4].replace(' ', '_')
#         genome_to_search_in = argv[5]
#
#         breakpoint_reverseComplement_contig = findPOI(strain_name_with_POI, genome_with_POI_file, anno_with_POI_file, product_of_interest,
#                              strain_name_to_search_in, genome_to_search_in)
#
#         if(len(argv) == 6):
#             ## no second annotation file given
#             rearrange_fasta_file(breakpoint_reverseComplement_contig[0], breakpoint_reverseComplement_contig[1], strain_name_to_search_in,
#                                  genome_to_search_in, product_of_interest, contig_name=breakpoint_reverseComplement_contig[2])
#
#         else:
#             ## rearrange the second annotation file
#             anno_to_search_in = argv[6]
#             rearrange_fasta_file(breakpoint_reverseComplement_contig[0], breakpoint_reverseComplement_contig[1], strain_name_to_search_in,
#                                  genome_to_search_in, product_of_interest, breakpoint_reverseComplement_contig[2])
#             rearrang_genbank_file(breakpoint_reverseComplement_contig[0], breakpoint_reverseComplement_contig[1], anno_to_search_in, strain_name_to_search_in)

# if(__name__ == '__main__'):
    # print (sys.argv)
    # main(sys.argv[1:])
    # para = ["6BC_Cpsittaci",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_genome.fasta",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_annotation_edit.gbk",
    #         "delta-aminolevulinic acid dehydratase"]

    # para = ["10DC88_Cavium",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_avium/10DC88_Cavium_genome.fasta",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_avium/10DC88_Cavium_annotation.gbk",
    #         "delta-aminolevulinic acid dehydratase",
    #         "11DC096_Cavium",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_avium/11DC096_Cavium_chromosome_node1.fasta"]

    # para = ["6BC_Cpsittaci",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_genome.fasta",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/6BC_Cpsittaci_annotation_edit.gbk",
    #         "delta-aminolevulinic acid dehydratase",
    #         "02DC15_Cpsittaci",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/02DC15_Cpsittaci_genome.fasta",
    #         "/mnt/prostlocal/marie/chlamydia_comparison/chlamydia_psittaci/02DC15_Cpsittaci_annotation.gbk"]
    #
    # main(para)
