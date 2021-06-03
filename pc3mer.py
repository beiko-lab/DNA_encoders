import numpy as np
from my_pseknc import Pseknc
from base_for_encoders import BaseForEncoder
import os


class Pc3mer(Pseknc, BaseForEncoder):
    def __init__(self,
                 classes=None,
                 max_samples=1008,
                 combine_class_samples=True,
                 folder_for_output=None):
        if classes:
            binary = [False if len(classes) != 2 else True for classes in [classes]][0]
        else:
            binary = False

        super().__init__(binary=binary)

        self.folder_for_output = folder_for_output
        self.path_db = os.path.join(self.folder_for_output, "Pc3mer")
        os.makedirs(self.path_db, exist_ok=True)
        self.combine_class_samples = combine_class_samples
        self.max_samples_per_class = max_samples
        self.binary = binary
        self.file_name_prefix = ""

        if classes:
            self.classes = np.array([self.synonyms(c) for c in classes])
            if binary:
                assert len(self.classes) == 2
        else:
            self.classes = None

    def property_values_along_the_sequence(self, seq, prop_name):
        """
        Finds all the kmers that compose the sequence keeping the order they appeared in the sequence
        Then, finds the property's values for each one of these k-mers
        :param seq: dna sequence to be encoded
        :param prop_name: name of the property. self.properties = ['Bendability-DNAse', 'Bendability-consensus', 'Trinucleotide GC Content', 'Nucleosome positioning',
         'Consensus_roll', 'Consensus_Rigid', 'Dnase I', 'Dnase I-Rigid', 'MW-Daltons', 'MW-kg', 'Nucleosome',
         'Nucleosome-Rigid']
        :return:
        """
        seq_into_oligonucs_index, _ = self.encode_into_kmers(seq)
        prop_as_feature = np.array(self.supInfo[prop_name]).astype(float)[
            seq_into_oligonucs_index
        ]
        return prop_as_feature

    def encode_fasta_into_pseknc_physicochemical_properties(
        self, input_path, prop_name, L=58, append_seq_id=True, verbose=False
    ):
        file_in = open(input_path, "r")
        sigma = ""
        feature_vector = []

        for i, line in enumerate(file_in):
            line_encoded = []
            if line == "":
                continue
            elif line[0] == ">":
                sigma = line[1:-1]
            else:
                if line[-1:] == "\n":
                    line58 = line[:-1]
                else:
                    line58 = line

                if len(line58) != L:
                    print(
                        "Line {} has size {} instead of {}.".format(i, len(line58), L)
                    )
                    continue

                line_encoded = list(
                    self.property_values_along_the_sequence(line58, prop_name)
                )
                if verbose:
                    print("line length {}".format(len(line58)))
                    print("encoded line length {}".format(len(line_encoded)))

                if append_seq_id:
                    line_encoded.append(self.synonyms(sigma))

            if line_encoded != []:
                feature_vector.append(line_encoded)

        file_in.close()

        return feature_vector
