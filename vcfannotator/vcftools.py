import collections
import vcf

class VCFTools(object):
    def __init__(self, vcf_file):
        self.reader = vcf.Reader(open(vcf_file, 'r'))
        self.records = [record for record in self.reader]
        
    def add_information_fields(self, info_list):
        for i in info_list:
            id, num, type, desc = i
            self.__add_information_field(id, num, type, desc)
            
    def __add_information_field(self, id, num, type, desc):
        _Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc'])
        if id:
            self.reader.infos[id] = _Info(id, num, type, desc)

    def write_vcf(self, output='/dev/stdout'):
        vcf_writer = vcf.Writer(open(output, 'w'), self.reader)
        for record in self.records:
            vcf_writer.write_record(record)
