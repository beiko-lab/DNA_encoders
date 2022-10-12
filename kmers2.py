import pandas as pd
import numpy as np
from itertools import product
from git.beiko_lab.DNA_encoders.base_for_encoders import BaseForEncoder
from Codes.Libraries.raw_data.fasta import parse_fasta


class Kmers2(BaseForEncoder):
    """
    Decompose a sequence into kmers and count them.
    This version 2 is almost twice faster than version 1 (i.e. Class Kmers)
    """
    def __init__(self, binary=True, list_of_ks=[4]):
        super().__init__(binary=binary)
        self.list_of_ks = list_of_ks
        self.oligonucs = {
            1: ["A", "C", "G", "T"],
            2: ["".join(c) for c in product("ACGT", repeat=2)],
            3: ["".join(c) for c in product("ACGT", repeat=3)],
            4: ["".join(c) for c in product("ACGT", repeat=4)],
            5: ["".join(c) for c in product("ACGT", repeat=5)],
            6: ["".join(c) for c in product("ACGT", repeat=6)],
        }

    def bulk_count_kmers_several_ks(self, seqs, list_of_ks=None):
        """
        Parameters
        ----------
        seq : str
            DNA sequence to encode
        list_of_ks : list
                List of values for the k-mers encoder. The values must be integers.

        Returns
        ----------
        feature_vector : np.ndarray
            Array of all the encoded sequences.
        """
        if list_of_ks is None:
            list_of_ks = self.list_of_ks

        # decompose each seq in series into the ks in list of ks
        decompose = []
        length = len(seqs[0])
        for seq in seqs:
            d=[]
            for k in list_of_ks:
                d.extend([[seq[i:i + k]] for i in range(0, length - k + 1)])
            # converting list to array and counting unique elements
            u,c = np.unique(d, return_counts=True)
            decompose.append(pd.DataFrame(c.reshape([1,u.shape[0]]), columns=u))

        del seqs

        # create an empty DataFrame
        columns = []
        for k in list_of_ks:
            columns.extend(self.oligonucs[k])

        feature_vector = pd.DataFrame(columns=columns)

        for d in decompose:
            feature_vector = pd.concat([feature_vector, d], ignore_index=True, axis=0, sort=True)

        feature_vector.fillna(value=0, inplace=True)

        return feature_vector

    # def count_kmers_single_k(self, seq, k=3):
    #     """
    #     Count k-mer occurrences in a given seq.
    #
    #     Parameters
    #     ----------
    #     seq : string
    #         A single DNA sequence.
    #     k : int
    #         The value of k for which to count kmers.
    #
    #     Returns
    #     -------
    #     counts : pd.DataFrame
    #         DataFrame where columns are oligonucleotides and line 0 has the # of times each oligonucleotide appears in seq
    #         A dictionary of counts keyed by their individual kmers (strings
    #         of length k).
    #     """
    #     # Calculate how many kmers of length k there are
    #     num_kmers = len(seq) - k + 1
    #     # decompose sequence into k-mers
    #     decomposed = [seq[i:i + k] for i in range(0, num_kmers)]
    #     # count how many k-mers there are in the list 'decomposed'
    #     u,c = np.unique(decomposed, return_counts=True)
    #     df = pd.DataFrame(c, columns=u)
    #     # Return the final counts
    #     return df

    def encode_fasta_file(
            self,
            fastafile: str,
            outputfile: str,
            list_of_ks: list = None,
            add_class_to_entries: bool = True,
    ):
        """
        Encodes a fastafile into a feature vector of k-mers for k in list_of_k
        and stores it into outputfile.

        Parameters
        ----------
        fastafile : string
            Each entry in this file must have 58 bp
        list_of_ks : list
            List of the numbers we want k to have
        outputfile : string
            Path and file name where to store the encoded sequences
        add_class_to_entries : boolean
            If the sequence's label should be appended to the encoded sequence
        verbose : boolean
            If the print statements should be enabled

        Returns
        -------
        None

        """
        if not list_of_ks:
            list_of_ks = self.list_of_ks
        else:
            self.list_of_ks = list_of_ks

        print(f"Encoding {fastafile} into k-mers for k in {list_of_ks}...")

        regs = parse_fasta(fastafile, split_key_into_keydescription=False, notify_if_more_than_one_contig=False)

        # count the number of oligonucleotides
        kmers_to_save_in_file = self.bulk_count_kmers_several_ks(regs.sequence.values)

        if add_class_to_entries:
            kmers_to_save_in_file = pd.concat([kmers_to_save_in_file, regs[["key"]]], axis=1)

        del regs
        kmers_to_save_in_file.to_csv(outputfile, header=False, index=False)
        print(f"Encoded sequence saved at {outputfile}")