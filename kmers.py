import pandas as pd
import numpy as np
from itertools import product

class Kmers():

    def __init__(self):
        self.oligonucs = {
            1: ["A","C","G","T"],
            2: [''.join(c) for c in product('ACGT', repeat=2)],
            3: [''.join(c) for c in product('ACGT', repeat=3)],
            4: [''.join(c) for c in product('ACGT', repeat=4)],
            5: [''.join(c) for c in product('ACGT', repeat=5)],
            6: [''.join(c) for c in product('ACGT', repeat=6)]
        }

    def several_k_mers(self,seq,ks):
        """

        :param seq:
        :param ks: is a list of k's
        :return:
        """
        feature_vector = np.array([])
        for k in ks:
            feature_vector = np.concatenate((feature_vector,self.count_kmers(seq, k=k).values),axis=None)
        return feature_vector.flatten()

    def count_kmers(self, seq, k=3):
        """Count kmer occurrences in a given seq.
        Parameters
        ----------
        seq : string
            A single DNA sequence.
        k : int
            The value of k for which to count kmers.

        Returns
        -------
        counts : DataFrame, columns are oligonucleotides and line 0 has the # of times each oligonucleotide appears in seq
            A dictionary of counts keyed by their individual kmers (strings
            of length k).
        """

        oligonuc = self.oligonucs[k]
        row = np.array([0] * len(oligonuc))
        row.shape = (1, len(oligonuc))
        dataFrame = pd.DataFrame(row, columns=oligonuc)

        # Calculate how many kmers of length k there are
        num_kmers = len(seq) - k + 1

        # Loop over the kmer start positions
        for i in range(num_kmers):
            # Slice the string to get the kmer
            kmer = seq[i:i + k]
            # Increment the count for this kmer
            dataFrame.loc[0, kmer] += 1

        # Return the final counts
        return dataFrame

    def convertFastaIntoKmers(self, fastafile, k, output):
        fastafile = open(fastafile, "r")
        kmerResFile = open(output, "w")
        sigma_class = None
        kmers_to_file = ''
        for line in fastafile:

            # create and prepare the dictionary
            if (line[0] == '>'):
                sigma_class = line[1:]
            else:
                # count the number of oligonucleotides
                kmersCount=self.count_kmers(line[:-1],k)
                # save the count into a file
                kmers_to_file +=str(kmersCount.values[0]).replace('  ',' ').replace('[ ','[').replace('[','').\
                                  replace('\n','').replace(' ',','). replace(',,',','). replace(']',",")+sigma_class+'\n'

        kmerResFile.write(kmers_to_file)
        fastafile.close()
        kmerResFile.close()


    def convertFastaIntoSeveralKmers(self, fastafile, list_of_k , outputfile, add_class_to_entries=True):
        """
        Encodes a fastafile into a feature vector of several k-mers for k in list_of_k.

        :param fastafile: Each entry in this file must have 58 bp
        :param list_of_k: List of the numbers we want k to have
        :param outputfile:
        :param add_class_to_entries:
        :return:
        """
        fastafile = open(fastafile ,"r")
        sigma=''
        kmers_to_save_in_file = pd.DataFrame([])
        with open(outputfile, 'w') as f:
            for line in fastafile:
                if (line[0] == '>'):
                    sigma = line[1:-1]
                else:
                    # count the number of oligonucleotides
                    kmersCount=self.several_k_mers(line[:-1],list_of_k)
                    if add_class_to_entries:
                        kmersCount = np.concatenate((kmersCount,np.asarray(sigma).flatten()))
                    # save the counts for all seqs in fastafile into a file
                    # but first keep it in memory till the entire file is encoded
                    a=pd.DataFrame(data=kmersCount).transpose()
                    kmers_to_save_in_file = pd.concat([kmers_to_save_in_file, a], axis=0, ignore_index=True)
            kmers_to_save_in_file.to_csv(f, header=False, index=False)
        fastafile.close()