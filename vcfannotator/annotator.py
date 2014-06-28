from vcfannotator import genbank
from vcfannotator import vcftools
from Bio.Seq import Seq

class Annotator(object):
    def __init__(self, gb_file=False, vcf_file=False):
        self.__gb = genbank.GenBank(gb_file)
        self.__vcf = vcftools.VCFTools(vcf_file)
        self.__gb.base_by_pos(1731)
        self.__gb.base_by_pos(1732)
       
    def annotate_vcf_records(self):
        for record in self.__vcf.records:
            self.__gb.index = record.POS
            record_annotation = self.__gb.annotation_by_position()
            
            if len(record.ALT) > 1:
                print 'Multiple SNPs ({0}->{1})\t{2}'.format(record.REF, record.ALT,record_annotation)
            elif len(record.ALT[0]) > 1 or len(record.REF) > 1:
                print 'Indel ({0}->{1})\t{2}'.format(record.REF, record.ALT[0],record_annotation)
            else:
                if record_annotation[1] is not None:
                    codon_info = self.__gb.codon_by_position(record.POS)
                    if record_annotation[0] == -1:
                        record.REF = Seq(str(record.REF)).complement()
                        record.ALT[0] = Seq(str(record.ALT[0])).complement()
                    alt_codon = list(codon_info[0])
                    alt_codon[codon_info[1]] = str(record.ALT[0])
                    ref_aa = codon_info[0].translate()
                    alt_aa = Seq(''.join(alt_codon)).translate()
                    
                    print '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}{7}{8}'.format(record.POS, record.REF, 
                                                 record.ALT, codon_info[0], ''.join(alt_codon), 
                                                 codon_info[1], ref_aa, codon_info[2], alt_aa)
                else:
                    print record.POS, record.REF, record.ALT, 'Intergenic'

