# DNA encoders

This package can encode DNA sequences into:

* Pc3mer
* Pc3mer stats
* PseKNC
* K-mers

The DNA sequence must be in [fasta](https://en.wikipedia.org/wiki/FASTA_format) format with lines of 58 nucleotides and extension _.fna_ or _.fasta_. Example:

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

After encoding the sequences, the algorithm stores them at _folder_for_output_/output/_encoder_name_/.

## Encoders
* [Pc3mer](#pc3mer)
* [Pc3mer stats](#pc3mer-stats)
* [PseKNC](#pseknc)
* [K-mers](#k-mers)

## Pc3mer
Using a table [[1]](#1) with twelve physicochemical properties values for each 3-mers, we standardize the values and calculate pc3mer by decomposing the input sequence into 3-mers and replacing each 3-mers in the order they appear in the sequence by its value for a given physicochemical property.

The physicochemical properties are Bendability-DNAse, Bendability-consensus, Trinucleotide GC Content, Nucleosome positioning, Consensus_roll, Consensus_Rigid, Dnase I, Dnase I-Rigid, MW-Daltons, MW-kg, Nucleosome, and Nucleosome-Rigid [[1]](#1).

Example 1:
* Sequence: GGGA...
* (Algorithm's step 1) It decomposes the sequence into 3-mers: GGG, GGA, ...
* Standardized physicochemical properties for GGG: 
  * 'Bendability-DNAse': 0.07230364 
  * 'Bendability-consensus': 0.3577835
  * 'Trinucleotide GC Content': 1.73205081
  * Etc.
* Standardized physicochemical properties for GGA:
  * 'Bendability-DNAse': 0.26511335
  * 'Bendability-consensus': -0.0969693
  * 'Trinucleotide GC Content': 0.57735027
  * Etc.
* (Algorithm's step 2) For each property, it replaces each 3-mer by its value for that property and stores as _Pc3mer/<property_name>.md_. 
  * For Bendability-DNAse, the encoded sequence will look like: ```[0.07230364, 0.26511335, ..., sample_class]```
  * For Bendability-consensus, the encoded sequence will look like: ```[0.3577835, -0.0969693, ..., sample_class]```
  * For Trinucleotide GC Content, the encoded sequence will look like: ```[1.73205081, 0.57735027, .., sample_class]```
  * Etc.
* When combining the individual files, make sure to delete the _sample_class_ column for all properties but the last one, so the encoded sequence look like:
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
import os

input_fasta = "input_file.fasta" # path + file name
folder_for_output = os.path.join(os.getcwd(), "output") # path to store the output

encoder = FastaToPc3mer(folder_for_output=folder_for_output)
encoder.convert_to_pc3mer(input_fasta)
```

Output:
It creates the folder "Pc3mer" in _<folder_for_output>_ and twelve _.md_ files, each one containing the encoded sequences for one of the properties. Examples: output/Pc3mer/Bendability-consensus.md and output/Pc3mer/Dnase I.md.

## Pc3mer stats
Encodes a sequence into pc3mer and then get a set of statistics over the encoded sequence. The statistics are:
* minimum, 
* maximum,
* mean,
* standard deviation,
* median, and 
* variance. 

Usage:
```python
from fasta_to_pc3mer import FastaToPc3mer
import os

input_fasta = "input_file.fasta" # path + file name
folder_for_output = os.path.join(os.getcwd(), "output") # path to store the output

encoder = FastaToPc3mer(folder_for_output=folder_for_output)
encoder.convert_to_pc3mer_stats(input_fasta)
```

## PseKNC
This implementation of PseKNC I [[1]](#1), [[2]](#2) decomposes the sequence into 3-mers and maps them to physicochemical property values specific for each word that is used to calculate scores. The scores, called Theta_{i}, are concatenated to the 3-mer decomposition and refers to all 3-mers _i_ nucleotides distant from each other, for _i_ in [1, 2]. The final array is devided by the sum of 3-mers counts and Theta scores.

Usage:
```python
from my_pseknc import Pseknc
import os

input_fasta = "input_file.fasta" # path + file name
folder_for_output = os.path.join(os.getcwd(), "output") # path to store the output
outputFile = os.path.join(folder_for_output, "/pseknc.csv")

encoder = Pseknc()
encoder.encode_fasta_into_pseknc(input_fasta, outputFile)
```


## k-mers
It counts the frequency of k-mers, an enumeration of all “words” of length _k_, for k in a given interval, in the DNA sequence. 

Example 1: 
* For _k_ in [2, 2] there are 4^k = 4^2 = 16 possible words: 
  ```{AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT}```
* Sequence: GGGA
* Decomposing into k-mers: GG, GG, GA
* Algorithm's output:
  ```[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0]```


Example 2:
* For _k_ in [1, 2] there are 4^1 + 4^2 = 4 + 16 possible words: 
  ```{A, C, G, T, AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT}``` 
* Sequence: GGGA
* Decomposed into k-mers: G, G, G, A, GG, GG, GA
* Algorithm's output:
  ```[1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0]```
  
Usage:
```python
from kmers import Kmers
import os

# defining list of ks
k_start = 1
k_end = 5
k_values = list(range(k_start, k_end+1))

input_fasta = "input_file.fasta" # path + file name
folder_for_output = os.path.join(os.getcwd(), "output") # path to store the output
outputFile = os.path.join(folder_for_output, f"/{k_start}_to_{k_end}_mers.csv")

encoder = Kmers()
encoded = encoder.convertFastaIntoSeveralKmers(input_fasta, k_values, outputFile)
```

## References
<a id="1">[1]</a> W. Chen, T. Y. Lei, D. C. Jin, H. Lin, and K. C. Chou, “PseKNC: A flexible web server for generating pseudo K-tuple nucleotide composition,” Anal. Biochem., vol. 456, no. 1, pp. 53–60, 2014, doi: 10.1016/j.ab.2014.04.001.

<a id="2">[2]</a> W. Chen, X. Zhang, J. Brooker, H. Lin, L. Zhang, and K.-C. Chou, “PseKNC-General: a cross-platform package for generating various modes of pseudo nucleotide compositions,” Bioinformatics, vol. 31, no. 1, pp. 119–120, Jan. 2015, doi: 10.1093/bioinformatics/btu602.

