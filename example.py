from fasta_to_pc3mer import FastaToPc3mer
from my_pseknc import Pseknc
from kmers import Kmers
import os

input_file_path = "input_file.fasta"
folder_for_output = os.path.join(os.getcwd(), "output")

"""
Converting fasta file into physicochemical properties 
"""
encoder = FastaToPc3mer(folder_for_output=folder_for_output)
encoder.convert_to_pc3mer(input_file_path)

"""
Converting fasta file into physicochemical statistics 
"""
encoder = FastaToPc3mer(folder_for_output=folder_for_output)
encoder.convert_to_pc3mer_stats(input_file_path)


"""
Converting fasta file into PseKNC - I 
"""
outputFile = folder_for_output + "/pseknc.csv"

encoder = Pseknc()
encoder.encode_fasta_into_pseknc(input_file_path, outputFile)


"""
Converting fasta file into Several k-mers 
"""
k_values = list(range(1,7))
outputFile = folder_for_output + "/1_to_6_mers.csv"

encoder = Kmers()
encoded = encoder.convertFastaIntoSeveralKmers(input_file_path,k_values,outputFile)