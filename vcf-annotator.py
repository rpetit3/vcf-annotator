#! /usr/bin/env python3
"""
vcf-annotator.py [-h] [--output STRING] [--version] VCF_FILE GENBANK_FILE

Annotate variants from a VCF file using the reference genome's GenBank file.

positional arguments:
  VCF_FILE         VCF file of Variants
  GENBANK_FILE     GenBank file of the reference genome.

optional arguments:
  -h, --help       show this help message and exit
  --output STRING  File to write VCF output to (Default STDOUT).
  --version        show program's version number and exit

Example Usage:
./vcf-annotator.py example-data/example.vcf example-data/example.gb
"""
import collections
from Bio import SeqIO
from Bio.Seq import Seq
import vcf
VERSION = 0.5


class Annotator(object):
    """Annotate a given VCF file according to the reference GenBank."""

    def __init__(self, gb_file=False, vcf_file=False):
        """Initialize variables."""
        self.__annotated_features = ["CDS", "tRNA", "rRNA", "ncRNA",
                                     "misc_feature"]
        self.__gb = GenBank(gb_file)
        self.__vcf = VCFTools(vcf_file)
        self.add_annotation_info()

    def add_annotation_info(self):
        """Add custom VCF info fields."""
        self.__vcf.add_information_fields([
            ['RefCodon', None, 'String', 'Reference codon'],
            ['AltCodon', None, 'String', 'Alternate codon'],
            ['RefAminoAcid', None, 'String', 'Reference amino acid'],
            ['AltAminoAcid', None, 'String', 'Alternate amino acid'],
            ['CodonPosition', '1', 'Integer', 'Codon position in the gene'],
            ['SNPCodonPosition', '1', 'Integer', 'SNP position in the codon'],
            ['AminoAcidChange', None, 'String', 'Amino acid change'],
            ['IsSynonymous', '1', 'Integer',
             '0:nonsynonymous, 1:synonymous, 9:N/A or Unknown'],
            ['IsTransition', '1', 'Integer',
             '0:transversion, 1:transition, 9:N/A or Unknown'],
            ['IsGenic', '1', 'Integer', '0:intergenic, 1:genic'],
            ['IsPseudo', '1', 'Integer', '0:not pseudo, 1:pseudo gene'],
            ['LocusTag', None, 'String', 'Locus tag associated with gene'],
            ['Gene', None, 'String', 'Name of gene'],
            ['Note', None, 'String', 'Note associated with gene'],
            ['Inference', None, 'String', 'Inference of feature.'],
            ['Product', None, 'String', 'Description of gene'],
            ['ProteinID', None, 'String', 'Protein ID of gene'],
            ['Comments', None, 'String', 'Example: Negative strand: T->C'],
            ['VariantType', None, 'String', 'Indel, SNP, Ambiguous_SNP'],
            ['FeatureType', None, 'String', 'The feature type of variant.'],
        ])

    def annotate_vcf_records(self):
        """Annotate each record in the VCF acording to the input GenBank."""
        for record in self.__vcf.records:
            self.__gb.index = record.POS

            # Set defaults
            record.INFO['RefCodon'] = '.'
            record.INFO['AltCodon'] = '.'
            record.INFO['RefAminoAcid'] = '.'
            record.INFO['AltAminoAcid'] = '.'
            record.INFO['CodonPosition'] = '.'
            record.INFO['SNPCodonPosition'] = '.'
            record.INFO['AminoAcidChange'] = '.'
            record.INFO['IsSynonymous'] = 9
            record.INFO['IsTransition'] = 9
            record.INFO['Comments'] = '.'
            record.INFO['IsGenic'] = '0'
            record.INFO['IsPseudo'] = '0'
            record.INFO['LocusTag'] = '.'
            record.INFO['Gene'] = '.'
            record.INFO['Note'] = '.'
            record.INFO['Inference'] = '.'
            record.INFO['Product'] = '.'
            record.INFO['ProteinID'] = '.'
            record.INFO['FeatureType'] = 'inter_genic'

            # Get annotation info
            if self.__gb.feature_exists:
                record.INFO['FeatureType'] = self.__gb.feature.type
                if self.__gb.feature.type in self.__annotated_features:
                    feature = self.__gb.feature
                    if feature.type == "CDS":
                        record.INFO['IsGenic'] = '1'

                    qualifiers = {
                        'Note': 'note', 'LocusTag': 'locus_tag',
                        'Gene': 'gene', 'Product': 'product',
                        'ProteinID': 'protein_id',
                        'Inference': 'inference'
                    }

                    if feature.type == "tRNA":
                        qualifiers['Note'] = 'anticodon'
                    for k, v in qualifiers.items():
                        if v in feature.qualifiers:
                            # Spell out semi-colons, commas and spaces
                            record.INFO[k] = feature.qualifiers[v][0].replace(
                                ';', '[semi-colon]'
                            ).replace(
                                ',', '[comma]'
                            ).replace(
                                ' ', '[space]'
                            )
                            if v == 'anticodon':
                                record.INFO[k] = 'anticodon{0}'.format(
                                    record.INFO[k]
                                )

                    if 'pseudo' in feature.qualifiers:
                        record.INFO['IsPseudo'] = '1'

            # Determine variant type
            if record.is_indel:
                if record.is_deletion:
                    record.INFO['VariantType'] = 'Deletion'
                else:
                    record.INFO['VariantType'] = 'Insertion'
            else:
                if len(record.ALT) > 1:
                    record.ALT = self.__gb.determine_iupac_base(record.ALT)
                    record.INFO['VariantType'] = 'Ambiguous_SNP'
                else:
                    if record.is_transition:
                        record.INFO['IsTransition'] = 1
                    else:
                        record.INFO['IsTransition'] = 0
                    record.INFO['VariantType'] = 'SNP'

                if int(record.INFO['IsGenic']):
                    alt_base = str(record.ALT[0])

                    # Determine codon information
                    codon = self.__gb.codon_by_position(record.POS)
                    record.INFO['RefCodon'] = ''.join(list(codon[0]))
                    record.INFO['SNPCodonPosition'] = codon[1]
                    record.INFO['CodonPosition'] = codon[2]

                    # Adjust for ambiguous base and negative strand.
                    if feature.strand == -1:
                        alt_base = str(
                            Seq(alt_base).complement()
                        )

                        record.INFO['Comments'] = 'Negative:{0}->{1}'.format(
                            Seq(record.REF).complement(),
                            alt_base
                        )

                    # Determine alternates
                    record.INFO['AltCodon'] = list(record.INFO['RefCodon'])
                    record.INFO['AltCodon'][
                        record.INFO['SNPCodonPosition']
                    ] = alt_base
                    record.INFO['AltCodon'] = ''.join(record.INFO['AltCodon'])
                    record.INFO['RefAminoAcid'] = Seq(
                        record.INFO['RefCodon']
                    ).translate()
                    record.INFO['AltAminoAcid'] = Seq(
                        record.INFO['AltCodon']
                    ).translate()
                    record.INFO['AminoAcidChange'] = '{0}{1}{2}'.format(
                        str(record.INFO['RefAminoAcid']),
                        record.INFO['CodonPosition'],
                        str(record.INFO['AltAminoAcid'])
                    )

                    if record.INFO['VariantType'] != 'Ambiguous_SNP':
                        ref = str(record.INFO['RefAminoAcid'])
                        alt = str(record.INFO['AltAminoAcid'])
                        if ref == alt:
                            record.INFO['IsSynonymous'] = 1
                        else:
                            record.INFO['IsSynonymous'] = 0

    def write_vcf(self, output='/dev/stdout'):
        """Write the VCF to the specified output."""
        self.__vcf.write_vcf(output)


class GenBank(object):
    """A class for parsing GenBank files."""

    def __init__(self, gb=False):
        """Inititalize variables."""
        self.__gb = SeqIO.read(open(gb, 'r'), 'genbank')
        self._index = None
        self.feature = None
        self.features = ["CDS", "rRNA", "tRNA", "ncRNA", "repeat_region",
                         "misc_feature"]
        self.build_position_index()
        self.gene_codons = {}

    @property
    def index(self):
        """Postion index for features."""
        return self._index

    @index.setter
    def index(self, value):
        self._index = self.__position_index[value - 1]
        self.__set_feature()

    def build_position_index(self):
        """Create a index for feature by position."""
        self.__position_index = [None] * len(self.__gb.seq)
        for i in range(len(self.__gb.features)):
            if self.__gb.features[i].type in self.features:
                start = int(self.__gb.features[i].location.start)
                end = int(self.__gb.features[i].location.end)
                self.__position_index[start:end] = [i] * (end - start)

    def __set_feature(self):
        if self._index is None:
            self.feature_exists = False
            self.feature = None
        else:
            self.feature_exists = True
            self.feature = self.__gb.features[self._index]

    def codon_by_position(self, pos):
        """Retreive the codon given a postion of a CDS feature."""
        if self._index not in self.gene_codons:
            self.split_into_codons()
        gene_position = self.position_in_gene(pos)
        codon_position = gene_position // 3
        return [self.gene_codons[self._index][codon_position],
                gene_position % 3,
                codon_position + 1]

    def split_into_codons(self):
        """Split the complete CDS feature in to a list of codons."""
        start = self.feature.location.start
        end = self.feature.location.end
        seq = ''.join(list(self.__gb.seq[start:end]))

        if self.feature.strand == -1:
            seq = Seq(seq).reverse_complement()

        self.gene_codons[self._index] = [
            seq[i:i + 3] for i in range(0, len(seq), 3)
        ]

    def position_in_gene(self, pos):
        """Return a codon postion in a gene."""
        if self.feature.strand == 1:
            return pos - self.feature.location.start - 1
        else:
            return self.feature.location.end - pos

    def base_by_pos(self, pos):
        """Print the base by position."""
        print(self.__gb.seq[pos - 1])

    def determine_iupac_base(self, bases):
        """
        Determine the IUPAC symbol for a list of nucleotides.

        Source: https://en.wikipedia.org/wiki/Nucleic_acid_notation
        List elements are in this order: [A,C,G,T]
        """
        if len(bases) > 1:
            iupac_notation = {
                'W': [True, False, False, True],
                'S': [False, True, True, False],
                'M': [True, True, False, False],
                'K': [False, False, True, True],
                'R': [True, False, True, False],
                'Y': [False, True, False, True],
                'B': [False, True, True, True],
                'D': [True, False, True, True],
                'H': [True, True, False, True],
                'V': [True, True, True, False],
                'N': [False, False, False, False]
            }

            base_condition = [base in bases for base in ['A', 'C', 'G', 'T']]
            for symbol, iupac_condition in iupac_notation.items():
                if iupac_condition == base_condition:
                    return symbol

    def is_transition(self, ref_base, alt_base):
        """
        Identify SNP as being a transition or not.

        1: Transition, 0:Transversion
        """
        substitution = ref_base + alt_base
        transition = ['AG', 'GA', 'CT', 'TC']

        if substitution in transition:
            return 1

        return 0


class VCFTools(object):
    """A class for parsing VCF formatted files."""

    def __init__(self, vcf_file):
        """Initialize variables."""
        self.reader = vcf.Reader(open(vcf_file, 'r'))
        self.records = [record for record in self.reader]

    def add_information_fields(self, info_list):
        """Add a given list of information fields to the VCF."""
        for i in info_list:
            id, num, type, desc = i
            self.__add_information_field(id, num, type, desc)

    def __add_information_field(self, id, num, type, desc):
        """Add a given information field to the VCF."""
        _Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc'])
        if id:
            self.reader.infos[id] = _Info(id, num, type, desc)

    def write_vcf(self, output='/dev/stdout'):
        """Write the VCF to a given output file."""
        vcf_writer = vcf.Writer(open(output, 'w'), self.reader)
        for record in self.records:
            vcf_writer.write_record(record)


if __name__ == '__main__':
    import argparse as ap
    import os
    import sys

    parser = ap.ArgumentParser(
        prog='vcf-annotator.py',
        conflict_handler='resolve',
        description=("Annotate variants from a VCF file using the reference "
                     "genome's GenBank file.")
    )
    parser.add_argument('vcf', metavar="VCF_FILE", type=str,
                        help='VCF file of SNPs')
    parser.add_argument('gb', metavar="GENBANK_FILE", type=str,
                        help='GenBank file of the reference genome.')
    parser.add_argument('--output', metavar="STRING", type=str,
                        default='/dev/stdout',
                        help='File to write VCF output to (Default STDOUT).')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {0}'.format(VERSION))

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    # Verify input exists
    if not os.path.isfile(args.gb):
        print('Unable to locate GenBank file: {0}'.format(args.gb))
        sys.exit(1)
    elif not os.path.isfile(args.vcf):
        print('Unable to locate VCF file: {0}'.format(args.vcf))
        sys.exit(1)

    annotator = Annotator(gb_file=args.gb, vcf_file=args.vcf)
    annotator.annotate_vcf_records()
    annotator.write_vcf(args.output)