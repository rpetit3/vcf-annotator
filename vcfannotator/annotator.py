from vcfannotator import genbank
from vcfannotator import vcftools
from Bio.Seq import Seq

class Annotator(object):
    def __init__(self, gb_file=False, vcf_file=False, length=15):
        self.length = length
        self.__gb = genbank.GenBank(gb_file)
        self.__vcf = vcftools.VCFTools(vcf_file)
        #self.__gb.base_by_pos(1731)
        #self.__gb.base_by_pos(1732)
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
            ['Substitution', '1', 'Integer', '0:synonymous, 1:nonsynonymous'],
            ['IsGenic', '1', 'Integer', '0:Intergenic, 1:Genic'],
            ['LocusTag', None, 'String', 'Locus tag associated with gene'],
            ['DBXref', None, 'String', 'Database ids associated with gene'],
            ['Gene', None, 'String', 'Name of gene'],
            ['Note', None, 'String', 'Note associated with gene'],
            ['Product', None, 'String', 'Description of gene'],
            ['ProteinID', None, 'String', 'Protein ID of gene'],
            ['Comments', None, 'String', 'Example: Negative strand: T->C'],
            ['FlankingRegion', None, 'String', "Sequence of flanking region"],
            ['VariantType', None, 'String', "Indel, SNP, Multiple SNPs"],
        ])
        
    def print_vcf():
        self.__vcf.print_metadata()
        self.__vcf.print_filter()
        self.__vcf.print_format()
        self.__vcf.print_info()
        self.__vcf.print_contigs()

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
            record.INFO['Substitution'] = '.'
            record.INFO['Comments'] = '.'
            record.INFO['IsGenic'] = '0'
            record.INFO['LocusTag'] = '.'
            record.INFO['DBXref'] = '.'
            record.INFO['Gene'] = '.'
            record.INFO['Note'] = '.'
            record.INFO['Product'] = '.'
            record.INFO['ProteinID'] = '.'
            record.INFO['FlankingRegion'] = self.__gb.get_flanking_region(self.length)
            
            # Get annotation info
            feature = self.__gb.feature
            if self.__gb.feature_exists():
                record.INFO['IsGenic'] = '1'
                record.INFO['LocusTag'] = feature.qualifiers['locus_tag']
                record.INFO['DBXref'] = feature.qualifiers['db_xref']
                record.INFO['Gene'] = feature.qualifiers['gene']
                record.INFO['Note'] = feature.qualifiers['note']
                record.INFO['Product'] = feature.qualifiers['product']
                record.INFO['ProteinID'] = feature.qualifiers['protein_id']

            # Determine variant type
            if len(record.ALT) > 1:
                record.INFO['VariantType'] = 'Multiple SNPs'
            elif len(record.ALT[0]) > 1 or len(record.REF) > 1:
                record.INFO['VariantType'] = 'Indel'
            else:
                record.INFO['VariantType'] = 'SNP'
                if int(record.INFO['IsGenic']):
                    codon_info = self.__gb.codon_by_position(record.POS)
                    if feature.strand == -1:
                        record.REF = Seq(str(record.REF)).complement()
                        record.ALT[0] = Seq(str(record.ALT[0])).complement()
                    alt_codon = list(codon_info[0])
                    alt_codon[codon_info[1]] = str(record.ALT[0])
                    ref_aa = codon_info[0].translate()
                    alt_aa = Seq(''.join(alt_codon)).translate()
                    
            
                    print '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}{7}{8}'.format(
                        record.POS, record.REF, record.ALT, codon_info[0], 
                        ''.join(alt_codon), codon_info[1], ref_aa, 
                        codon_info[2], alt_aa
                    )
