import vcf

class VCFTools(object):
    def __init__(self, vcf_file):
        self.__vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        self.records = [record for record in self.__vcf_reader]
        
    def print_vcf_records(self):
        for record in self.records:
            print record
            
    def variant_positions(self):
        return [record.POS for record in self.records]