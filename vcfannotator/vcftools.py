import collections
import vcf

class VCFTools(object):
    def __init__(self, vcf_file):
        self.reader = vcf.Reader(open(vcf_file, 'r'))
        self.records = [record for record in self.reader]
        self.number_conversion = {
            '.': None,  # Unknown number of values
            'A': -1,  # Equal to the number of alleles in a given record
            'G': -2,  # Equal to the number of genotypes in a given record
        }
        self._Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc'])
        
    def print_vcf_records(self):
        for record in self.records:
            print record
            
    def variant_positions(self):
        return [record.POS for record in self.records]
        
        
    def add_information_fields(self, info_list):
        for i in info_list:
            id, num, type, desc = i
            self.__add_information_field(id, num, type, desc)
            
    def __add_information_field(self, id, num, type, desc):
        if id:
            self.reader.infos[id] = self._Info(id, num, type, desc)
        
    def convert_number(self, number):
        '''
            Covert PyVCF's way of storing Number fields back to VCF standard

            https://github.com/jamescasbon/PyVCF/blob/master/vcf/parser.py            
        '''
        if number is None:
            return '.'
        elif number == -1:
            return 'A'
        elif number == -2:
            return 'G'
        else:
            return number
          
    def print_metadata(self):
        for k,v in self.reader.metadata.items():
            print '##{0}={1}'.format(k,v)
            
    def print_filter(self):
        '''
        Example: ##FILTER=<ID=LowConf,Description="DP < 10|| QD < 2.5">
        '''
        for i in self.reader.filters:
            print '##FILTER=<ID={0},Description="{1}">'.format(
                self.reader.filters[i].id, self.reader.filters[i].desc
            )
            
    def print_format(self):
        '''
        Example: ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        '''
        for i in self.reader.formats:
            print '##FORMAT=<ID={0},Number={1},Type={2},Description="{3}">'.format(
                self.reader.formats[i].id, 
                self.convert_number(self.reader.formats[i].num), 
                self.reader.formats[i].type, 
                self.reader.formats[i].desc
            )
    def print_info(self):
        '''
        Example: ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
        '''
        for i in self.reader.infos:
            print '##INFO=<ID={0},Number={1},Type={2},Description="{3}">'.format(
                self.reader.infos[i].id, 
                self.convert_number(self.reader.infos[i].num), 
                self.reader.infos[i].type, 
                self.reader.infos[i].desc
            )
            
    def print_contigs(self):
        '''
        Example: ##contig=<ID=gi|289163350|ref|NC_013774.1|,length=1767>
        '''
        for i in self.reader.contigs:
            print '##contig=<ID={0},length={1}>'.format(
                self.reader.contigs[i].id, 
                self.reader.contigs[i].length 
            )
