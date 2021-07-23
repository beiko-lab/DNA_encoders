import pandas as pd
import numpy as np
from itertools import product
from beiko_lab.DNA_encoders.base_for_encoders import BaseForEncoder


class Kmers(BaseForEncoder):
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

    def count_kmers_several_ks(self, seq, list_of_ks=None):
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
        feature_vector = np.array([])

        if list_of_ks is None:
            list_of_ks = self.list_of_ks

        for k in list_of_ks:
            feature_vector = np.concatenate(
                (feature_vector, self.count_kmers_single_k(seq, k=k).values), axis=None
            )
        return feature_vector.flatten()

    def count_kmers_single_k(self, seq, k=3):
        """
        Count k-mer occurrences in a given seq.

        Parameters
        ----------
        seq : string
            A single DNA sequence.
        k : int
            The value of k for which to count kmers.

        Returns
        -------
        counts : pd.DataFrame
            DataFrame where columns are oligonucleotides and line 0 has the # of times each oligonucleotide appears in seq
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
            kmer = seq[i: i + k]
            # Increment the count for this kmer
            dataFrame.loc[0, kmer] += 1

        # Return the final counts
        return dataFrame

    def encode_fasta_file(
            self,
            fastafile: str,
            outputfile: str,
            list_of_ks: list = None,
            add_class_to_entries: bool = True,
            verbose: bool = False,
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

        fastafile = open(fastafile, "r")
        sigma = ""
        kmers_to_save_in_file = pd.DataFrame([])
        for line in fastafile:
            if line.startswith(">"):
                sigma = line[1:-1]
            else:
                # count the number of oligonucleotides
                kmersCount = self.count_kmers_several_ks(line[:-1])
                if add_class_to_entries:
                    kmersCount = np.concatenate(
                        (kmersCount, np.asarray(self.synonyms(sigma)).flatten())
                    )

                if verbose:
                    print(f"Seq: {line[:-1]}\nSeq_encoded: {kmersCount}")

                # save the counts for all seqs in fastafile into a file
                # but first keep it in memory until the entire file is encoded
                a = pd.DataFrame(data=kmersCount).transpose()
                kmers_to_save_in_file = pd.concat(
                    [kmers_to_save_in_file, a], axis=0, ignore_index=True
                )
        kmers_to_save_in_file.to_csv(outputfile, header=False, index=False)
        print(f"Encoded sequence saved at {outputfile}")
        fastafile.close()
