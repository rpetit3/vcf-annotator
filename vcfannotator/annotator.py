from vcfannotator import genbank
from vcfannotator import vcftools
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
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
                    'Note':'note', 'LocusTag':'locus_tag', 'DBXref':'db_xref',
                    'Gene':'gene', 'Product':'product', 'ProteinID':'protein_id'
                }
                for k,v in qualifiers.items():
                    if v in feature.qualifiers:
                        record.INFO[k] = feature.qualifiers[v]

            # Determine variant type
            #if len(record.ALT) > 1:
            #    self.__gb.determine_iupac_base(record.ALT)
            #    record.INFO['VariantType'] = 'Multiple_SNPs'
            if record.ALT[0] is None:
                record.INFO['VariantType'] = 'Indel'
            elif len(record.ALT[0]) > 1 or len(record.REF) > 1:
                record.INFO['VariantType'] = 'Indel'
            else:
                if len(record.ALT) > 1:
                    record.ALT = self.__gb.determine_iupac_base(record.ALT)
                    record.INFO['VariantType'] = 'Ambiguous_SNP'
                else:
                    record.INFO['IsTransition'] = self.__gb.is_transition(
                        str(record.REF), str(record.ALT[0])
                    )
                    record.INFO['VariantType'] = 'SNP'
                    
                if int(record.INFO['IsGenic']):
                    
                    alt_base = str(record.ALT[0])
                    ref_base = str(record.REF)

                    #Determine codon information
                    codon = self.__gb.codon_by_position(record.POS)
                    record.INFO['RefCodon'] = codon[0]
                    record.INFO['SNPCodonPosition'] = codon[1]
                    record.INFO['CodonPosition'] = codon[2]
                    
                    # Adjust for ambiguous base and negative strand.
                    
                    
                    if feature.strand == -1:
                        alt_base = str(
                            Seq(alt_base, IUPAC.ambiguous_dna).complement()
                        )

                        record.INFO['Comments'] = 'Negative {0} -> {1}'.format(
                            Seq(record.REF).complement(), 
                            alt_base
                        )

                    # Determine alternates
                    record.INFO['AltCodon'] = list(record.INFO['RefCodon'])
                    record.INFO['AltCodon'][record.INFO['SNPCodonPosition']] = alt_base
                    record.INFO['AltCodon'] = ''.join(record.INFO['AltCodon'])
                    record.INFO['RefAminoAcid'] = record.INFO['RefCodon'].translate()
                    record.INFO['AltAminoAcid'] = Seq(
                        record.INFO['AltCodon'], 
                        IUPAC.ambiguous_dna
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