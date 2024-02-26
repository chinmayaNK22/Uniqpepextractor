# Uniqpepextractor
Tryptic peptides unique to a protein/species from a list of multiple proteins/species fasta files can be extracted.

```
>python Uniqpepextractor.py path_to_proteomeDB/proteomes output_SSPPs/ 2 6 30 trypsin -AA C,M

SSPE: Species Specific Peptide Extractor. Peptides list from the species specific input fasta file (Proteome) will be
generated and only those peptide sequences specific to a input fasta file (i.e. Species proteome) will be considered
when compared with the other input proteome databases

positional arguments:
  -ip                   Path to the folder in which proteome fasta files are stored
  -op                   Path to the folder in which all result files will be stored
  -mc                   Maximum missed cleavage to be considered for in-silico digestion by the protease
  -ml                   Minimum length of the peptide do be considered
  -ML                   Maximum length of peptide to be considered
  -enz                  Protease enzyme for in-silico digestion (Only Trypsin, LysC and Chymotrypsin can be applied)

options:
  -h, --help            show this help message and exit
  -AA PEPTIDES_WITH_AA, --peptides_with_AA PEPTIDES_WITH_AA
                        Peptides with these amino acids will not be considered for the exctraction (e.x. C,M)
```
