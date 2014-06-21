import vcf

class VCFTools(object):
    def __init__(self, vcf_file):
        self.__vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    def print_vcf_records(self):
        for record in self.__vcf_reader:
            print record