from vcfannotator import genbank
from vcfannotator import vcftools
from Bio.Seq import Seq
import vcf


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
            ['Substitution', '1', 'Integer', '0:synonymous, 1:nonsynonymous'],
            ['IsGenic', '1', 'Integer', '0:intergenic, 1:genic'],
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
            record.INFO['FlankingRegion'] = '.'

            # Get annotation info
            if self.__gb.feature_exists:
                feature = self.__gb.feature
                record.INFO['IsGenic'] = '1'
                qualifiers = {
                    'Note':'note', 'LocusTag':'locus_tag', 'DBXref':'db_xref',
                    'Gene':'gene', 'Product':'product', 'ProteinID':'protein_id'
                }
                for k,v in qualifiers.items():
                    if v in feature.qualifiers:
                        record.INFO[k] = feature.qualifiers[v]

            # Determine variant type
            if len(record.ALT) > 1:
                record.INFO['VariantType'] = 'Multiple_SNPs'
            elif len(record.ALT[0]) > 1 or len(record.REF) > 1:
                record.INFO['VariantType'] = 'Indel'
            else:
                record.INFO['VariantType'] = 'SNP'
                if int(record.INFO['IsGenic']):
                    alt_base = str(record.ALT[0])
                    ref_base = str(record.REF)
                    record.INFO['FlankingRegion'] = self.__gb.get_flanking_region(
                        record.ALT[0], record.POS, self.length
                    )
                    
                    #Determine codon information
                    codon = self.__gb.codon_by_position(record.POS)
                    record.INFO['RefCodon'] = codon[0]
                    record.INFO['SNPCodonPosition'] = codon[1]
                    record.INFO['CodonPosition'] = codon[2]
                    
                    # Adjust for negative strand
                    if feature.strand == -1:
                        alt_base = str(Seq(str(record.ALT[0])).complement())
                        record.INFO['FlankingRegion'] = Seq(record.INFO['FlankingRegion']).reverse_complement()
                        record.INFO['Comments'] = 'Negative {0} -> {1}'.format(
                            Seq(record.REF).complement(), 
                            Seq(str(record.ALT[0])).complement()
                        )
                        
                    # Determine alternates
                    record.INFO['AltCodon'] = list(record.INFO['RefCodon'])
                    record.INFO['AltCodon'][record.INFO['SNPCodonPosition']] = alt_base
                    record.INFO['AltCodon'] = ''.join(record.INFO['AltCodon'])
                    record.INFO['RefAminoAcid'] = record.INFO['RefCodon'].translate()
                    record.INFO['AltAminoAcid'] = Seq(record.INFO['AltCodon']).translate()
                    record.INFO['AminoAcidChange'] = '{0}{1}{2}'.format(
                        str(record.INFO['RefAminoAcid']),
                        record.INFO['CodonPosition'],
                        str(record.INFO['AltAminoAcid'])
                    )
                    
                    if str(record.INFO['RefAminoAcid']) == str(record.INFO['AltAminoAcid']):
                        record.INFO['Substitution'] = 'SYN'
                    else:
                        record.INFO['Substitution'] = 'nSYN'
                        
    def write_vcf(self):
        self.__vcf.write_vcf()