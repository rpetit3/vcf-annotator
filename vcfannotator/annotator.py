from vcfannotator import genbank
from vcfannotator import vcftools

class Annotator(object):
    def __init__(self, gb_file=False, vcf_file=False):
        self.__gb = genbank.GenBank(gb_file)
        self.__gb.print_genbank_stats()
        self.__vcf = vcftools.VCFTools(vcf_file)
        self.__vcf.print_vcf_records()
