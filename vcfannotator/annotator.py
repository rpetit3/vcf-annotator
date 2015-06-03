from vcfannotator import genbank
from vcfannotator import vcftools
from Bio.Seq import Seq


class Annotator(object):
    def __init__(self, gb_file=False, vcf_file=False, length=15):
        self.length = length
        self.__gb = genbank.GenBank(gb_file)
        self.__vcf = vcftools.VCFTools(vcf_file)
        self.add_annotation_info()

    def add_annotation_info(self):
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
            ['LocusTag', None, 'String', 'Locus tag associated with gene'],
            ['DBXref', None, 'String', 'Database ids associated with gene'],
            ['Gene', None, 'String', 'Name of gene'],
            ['Note', None, 'String', 'Note associated with gene'],
            ['Product', None, 'String', 'Description of gene'],
            ['ProteinID', None, 'String', 'Protein ID of gene'],
            ['Comments', None, 'String', 'Example: Negative strand: T->C'],
            ['VariantType', None, 'String', 'Indel, SNP, Ambiguous_SNP'],
        ])

    def annotate_vcf_records(self):
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
            record.INFO['LocusTag'] = '.'
            record.INFO['DBXref'] = '.'
            record.INFO['Gene'] = '.'
            record.INFO['Note'] = '.'
            record.INFO['Product'] = '.'
            record.INFO['ProteinID'] = '.'

            # Get annotation info
            if self.__gb.feature_exists:
                feature = self.__gb.feature
                record.INFO['IsGenic'] = '1'
                qualifiers = {
                    'Note': 'note', 'LocusTag': 'locus_tag',
                    'DBXref': 'db_xref', 'Gene': 'gene', 'Product':
                    'product', 'ProteinID': 'protein_id'
                }
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
                    record.INFO['IsTransition'] = 1 if record.is_transition else 0
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
                    record.INFO['AltCodon'][record.INFO['SNPCodonPosition']] = alt_base
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

    def write_vcf(self):
        self.__vcf.write_vcf()
