from Bio import SeqIO

class GenBank(object):
    def __init__(self, gb=False):
        self.__gb = SeqIO.read(open(gb, 'r'), 'genbank')
        
    def print_genbank_stats(self):
        print "Name %s, %i features" % (self.__gb.name, len(self.__gb.features))
