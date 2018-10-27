[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/vcf-annotator/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/vcf-annotator/badges/downloads.svg)](https://anaconda.org/bioconda/vcf-annotator)
<!-- [![Docker Repository on Quay.io](https://quay.io/repository/biocontainers/vcf-annotator/status "Docker Repository on Quay.io")](https://quay.io/repository/biocontainers/vcf-annotator) -->

*vcf-annotator* uses the reference GenBank file to add more details to the variant calls in a VCF.

# vcf-annotator
Using a reference GenBank file, *vcf-annotator* adds biological annotations to variants in a VCF file. A full list of annotations is descibed below, but these include amino acid changes, gene information, synonymous vs nonsynonymous, locus tag information, among many more. 

### Added Annotations
For each mutation, if applicable, the following annotations are added to the INFO column of the VCF.

| Annotation | Description |
|------------|-------------|
| RefCodon | Reference codon |
| AltCodon | Alternate codon |
| RefAminoAcid | Reference amino acid |
| AltAminoAcid | Alternate amino acid |
| CodonPosition | Codon position in the gene |
| SNPCodonPosition | SNP position in the codon |
| AminoAcidChange | Amino acid change |
| IsSynonymous | 0:nonsynonymous, 1:synonymous, 9:N/A or Unknown |
| IsTransition | 0:transversion, 1:transition, 9:N/A or Unknown |
| IsGenic | 0:intergenic, 1:genic |
| IsPseudo | 0:not pseudo, 1:pseudo gene |
| LocusTag | Locus tag associated with gene |
| Gene | Name of gene |
| Note | Note associated with gene |
| Inference | Inference of feature. |
| Product | Description of gene |
| ProteinID | Protein ID of gene |
| Comments | Example: Negative strand: T->C |
| VariantType | Indel, SNP, Ambiguous_SNP |
| FeatureType | The feature type of variant. |

# Installation
### Requirements
* Python >= 3.4
* [BioPython](http://biopython.org/) >= 1.65 
* [PyVCF](https://github.com/jamescasbon/PyVCF) == 0.6.8

### Bioconda
*vcf-annotator* is available from BioConda
```
conda install -c bioconda vcf-annotator
```

### From Source
```
git@github.com:rpetit3/vcf-annottor.git
cd vcf-annottor
pip3 install -r requirements.txt
python3 vcf-annottor.py YOUR_VCF.vcf REFERENCE.gb
```

Nothing much else to it, just a simple to read in a VCF and GenBank file and output an annotated VCF. Feel free to drop it in your $PATH somewhere!

# Usage
*vcf-annotator* requires an uncompressed VCF file and the corresponding reference GenBank file. It then outputs the annotated variants, by default to STDOUT, but this can be changed on runtime.

### Usage Output
```
python3 vcf-annotator.py
usage: vcf-annotator.py [-h] [--output STRING] [--version]
                        VCF_FILE GENBANK_FILE

Annotate variants from a VCF file using the reference genome's GenBank file.

positional arguments:
  VCF_FILE         VCF file of variants
  GENBANK_FILE     GenBank file of the reference genome.

optional arguments:
  -h, --help       show this help message and exit
  --output STRING  File to write VCF output to (Default STDOUT).
  --version        show program's version number and exit
```

##### --version Output
```
python3 vcf-annotator.py --version
vcf-annotator.py 0.5
```

### Example Usage
A VCF and GenBank file are included in the *example-data* directory. You can use these two files to verify the script is working properly.
```
python3 vcf-annotator.py example-data/example.vcf example-data/example.gb
```

### Disclaimer
This script has been developed only for microbial variant analysis. I've only tested on VCF files output from GATK, but I would assume if the VCF format is followed other VCF files should work as well. Currently for a ~3mb genome with ~20k mutations it takes about 10s to annotate the VCF file. Based on this information, I'm not sure how well it would work on larger genomes (if it would even work at all!).
  
