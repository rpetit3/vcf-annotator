"""Parse a reference GenBank file."""
from Bio import SeqIO
from Bio.Seq import Seq


class GenBank(object):
    """A class for parsing GenBank files."""

    def __init__(self, gb=False):
        """Inititalize variables."""
        self.__gb = SeqIO.read(open(gb, 'r'), 'genbank')
        self._index = None
        self.feature = None
        self.features = ["CDS", "rRNA", "tRNA", "ncRNA", "repeat_region",
                         "misc_feature"]
        self.build_position_index()
        self.gene_codons = {}

    @property
    def index(self):
        """Postion index for features."""
        return self._index

    @index.setter
    def index(self, value):
        self._index = self.__position_index[value - 1]
        self.__set_feature()

    def build_position_index(self):
        """Create a index for feature by position."""
        self.__position_index = [None] * len(self.__gb.seq)
        for i in xrange(len(self.__gb.features)):
            if self.__gb.features[i].type in self.features:
                start = int(self.__gb.features[i].location.start)
                end = int(self.__gb.features[i].location.end)
                self.__position_index[start:end] = [i] * (end - start)

    def __set_feature(self):
        if self._index is None:
            self.feature_exists = False
            self.feature = None
        else:
            self.feature_exists = True
            self.feature = self.__gb.features[self._index]

    def codon_by_position(self, pos):
        """Retreive the codon given a postion of a CDS feature."""
        if self._index not in self.gene_codons:
            self.split_into_codons()
        gene_position = self.position_in_gene(pos)
        codon_position = gene_position / 3
        return [self.gene_codons[self._index][codon_position],
                gene_position % 3,
                codon_position + 1]

    def split_into_codons(self):
        """Split the complete CDS feature in to a list of codons."""
        start = self.feature.location.start
        end = self.feature.location.end
        seq = ''.join(list(self.__gb.seq[start:end]))

        if self.feature.strand == -1:
            seq = Seq(seq).reverse_complement()

        self.gene_codons[self._index] = [
            seq[i:i + 3] for i in range(0, len(seq), 3)
        ]

    def position_in_gene(self, pos):
        """Return a codon postion in a gene."""
        if self.feature.strand == 1:
            return pos - self.feature.location.start - 1
        else:
            return self.feature.location.end - pos

    def base_by_pos(self, pos):
        """Print the base by position."""
        print(self.__gb.seq[pos - 1])

    def determine_iupac_base(self, bases):
        """
        Determine the IUPAC symbol for a list of nucleotides.

        Source: https://en.wikipedia.org/wiki/Nucleic_acid_notation
        List elements are in this order: [A,C,G,T]
        """
        if len(bases) > 1:
            iupac_notation = {
                'W': [True, False, False, True],
                'S': [False, True, True, False],
                'M': [True, True, False, False],
                'K': [False, False, True, True],
                'R': [True, False, True, False],
                'Y': [False, True, False, True],
                'B': [False, True, True, True],
                'D': [True, False, True, True],
                'H': [True, True, False, True],
                'V': [True, True, True, False],
                'N': [False, False, False, False]
            }

            base_condition = [base in bases for base in ['A', 'C', 'G', 'T']]
            for symbol, iupac_condition in iupac_notation.items():
                if iupac_condition == base_condition:
                    return symbol

    def is_transition(self, ref_base, alt_base):
        """
        Identify SNP as being a transition or not.

        1: Transition, 0:Transversion
        """
        substitution = ref_base + alt_base
        transition = ['AG', 'GA', 'CT', 'TC']

        if substitution in transition:
            return 1
        else:
            return 0
