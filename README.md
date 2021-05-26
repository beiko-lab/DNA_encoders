# DNA encoders

This package has the following encoders for DNA sequences:

***

The DNA sequence must be in [fasta](https://en.wikipedia.org/wiki/FASTA_format) with lines of fixed length and extension 
_.fna_ or _.fasta_. Example:

```bash
>Sigma++
GGTTTATTGCCTTGCAGCTGGCGAGAGACGGTATTGCTCATGCACAAGCCTTGTTCAG
>Sigma24
TGCCCTGACTTCACCCCGCTGTGTCTGCTTTTCCCGACTATTCTTAATGAGCTTCGAT
>Sigma--
AATGTGGATAGATATGAATTATTTTTCTCCTTAAGGATCATCCGTTATTTGGGTCGTT
>Sigma70
CAGTTATTTACCTTACTTTACGCGCGCGTAACTCTGGCAACATCACTACAGGATAGCG
>Sigma++
AAAAAGTTATACGCGGTGGAAACATTGCCCGGATAGTCTATAGTCACTAAGCATTAAA
...
```

After choosing the wanted encoder, the files with encoded sequences will be stored at 
_folder_for_output_/output/_encoder_name_/

## Encoders
* [Pc3mer](#pc3mer)
* [Pc3mer stats](#pc3mer)
* [PseKNC](#pseknc)
* [k-mers](#k-mers)

## Pc3mer
Using (add ref)'s table mapping 3-mers to a value for each one of the twelve physicochemical properties, we 
standardize the values and calculate pc3mer by decomposing the input sequence into 3-mers and then replacing each 
3-mers by its value for a given physicochemical property.

The properties names are Bendability-DNAse, Bendability-consensus, Trinucleotide GC Content, Nucleosome positioning, 
Consensus_roll, Consensus_Rigid, Dnase I, Dnase I-Rigid, MW-Daltons, MW-kg, Nucleosome, Nucleosome-Rigid (add ref).

Example 1:
* Sequence: GGGA...
* (Step 1) Decomposing it into 3-mers: GGG, GGA, ...
* Standardized physicochemical properties for GGG: 
  * 'Bendability-DNAse': 0.07230364, 
  * 'Bendability-consensus': 0.3577835
  * 'Trinucleotide GC Content': 1.73205081
  * Etc.
* Standardized physicochemical properties for GGA:
  * 'Bendability-DNAse': 0.26511335, 
  * 'Bendability-consensus': -0.0969693
  * 'Trinucleotide GC Content': 0.57735027
  * Etc.
* (Step 2) For each property, replace each 3-mer by its value for that property and store as _Pc3mer/property.md_. 
  * For Bendability-DNAse, the encoded sequence starts with: ```[0.07230364, 0.26511335, ..., sample_class]```
  * For Bendability-consensus, the encoded sequence starts with: ```[0.3577835, -0.0969693, ..., sample_class]```
  * For Trinucleotide GC Content, the encoded sequence starts with: ```[1.73205081, 0.57735027, .., sample_class]```
  * Etc.
* When the individual files are combined, make sure to delete the _sample_class_ for all the properties, but
  the last one, so the encoded sequence look like:
  ```
   [0.07230364, 0.26511335, ..., 0.3577835, -0.0969693, ..., 1.73205081, 0.57735027, ..., sample_class]
  ```
  
Example 2:
* Entry in the fasta file: 
```bash
>Sigma++
GCTGAAAATACGTTGAACGCTTACCGTCGCGATCTGTCAATGATGGTGGAGTGGTTGC
```

* Sequence encoded with Bendability-consensus:
```
0.919537039821516,0.919537039821516,1.3475397297384397,-0.605222559882525,-2.745236029467144,-2.745236029467144, -2.584735019498298,0.5717848498890153,-0.0702191899863703,0.06353164998766837,0.06353164998766837,..., Sigma++
```

Usage:
```python
from fasta_to_pc3mer import FastaToPc3mer
from my_pseknc import Pseknc
from kmers import Kmers
import os

input_fasta = "input_file.fasta" # path + file name
folder_for_output = os.path.join(os.getcwd(), "output") # path to store the output
encoder = FastaToPc3mer(folder_for_output=folder_for_output)
encoder.convert_to_pc3mer(input_fasta)
```

Output:
It creates the folder "Pc3mer" in _folder_for_output_ and twelve _.md_ files, each one containing the fasta file encoded
for one of the properties. Examples: output/Pc3mer/Bendability-consensus.md and output/Pc3mer/Dnase I.md.


## k-mers




