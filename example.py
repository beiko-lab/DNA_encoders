from fasta_to_pc3mer import FastaToPc3mer
from my_pseknc import Pseknc
from kmers import Kmers
import os

input_fasta = "input_file.fasta"
folder_for_output = os.path.join(os.getcwd(), "output")

"""
Converting fasta file into pc3mer properties 
"""
print("*pc3mer*")
encoder = FastaToPc3mer(folder_for_output=folder_for_output)
encoder.convert_to_pc3mer(input_fasta)

"""
Converting fasta file into pc3mers statistics 
"""
print("*pc3mer stats*")
encoder = FastaToPc3mer(folder_for_output=folder_for_output)
encoder.convert_to_pc3mer_stats(input_fasta)


"""
Converting fasta file into PseKNC - I 
"""
print("*PseKNC I*")
outputFile = os.path.join(folder_for_output, "pseknc.csv")
encoder = Pseknc()
encoder.encode_fasta_into_pseknc(input_fasta, outputFile)


"""
Converting fasta file into k-mers 
"""
print("*k-mers*")
# defining list of ks
k_start = 1
k_end = 2
k_values = list(range(k_start, k_end + 1))

outputFile = os.path.join(folder_for_output, f"{k_start}_to_{k_end}_mers.csv")

encoder = Kmers()
encoded = encoder.convertFastaIntoSeveralKmers(input_fasta, k_values, outputFile)
