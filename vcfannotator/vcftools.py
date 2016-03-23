"""Parse VCF files."""
import collections
import vcf


class VCFTools(object):
    """A class for parsing VCF formatted files."""

    def __init__(self, vcf_file):
        """Initialize variables."""
        self.reader = vcf.Reader(open(vcf_file, 'r'))
        self.records = [record for record in self.reader]

    def add_information_fields(self, info_list):
        """Add a given list of information fields to the VCF."""
        for i in info_list:
            id, num, type, desc = i
            self.__add_information_field(id, num, type, desc)

    def __add_information_field(self, id, num, type, desc):
        """Add a given information field to the VCF."""
        _Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc'])
        if id:
            self.reader.infos[id] = _Info(id, num, type, desc)

    def write_vcf(self, output='/dev/stdout'):
        """Write the VCF to a given output file."""
        vcf_writer = vcf.Writer(open(output, 'w'), self.reader)
        for record in self.records:
            vcf_writer.write_record(record)
