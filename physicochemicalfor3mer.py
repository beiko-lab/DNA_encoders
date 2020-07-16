import pandas as pd
import numpy as np
from my_pseknc import Pseknc
from base_for_encoders import BaseForEncoder
import os.path as osp
import os

class PhysicochemicalFor3mer(Pseknc, BaseForEncoder):


    def __init__(self, classes=None, max_samples=1008, combine_class_samples=True,
                 binary=True, folder_for_output=None):
        super().__init__()
        self.folder_for_output = folder_for_output
        self.path_db = self.folder_for_output + "Physicochemical_props/"
        if not osp.isdir(self.path_db):
            os.mkdir(self.path_db)
        self.combine_class_samples = combine_class_samples

        if classes:
            self.classes = np.array([self.synonimous(c) for c in classes])
        else:
            self.classes = None

        self.max_samples_per_class = max_samples
        self.binary = binary
        if binary and classes:
            assert len(self.classes) == 2

        self.file_name_prefix = ""

    def read_from_file(self):

        promfeatures = None

        for prop in self.properties:
            features = pd.io.parsers.read_csv(self.path_db + self.file_name_prefix + prop + ".md", header=None, \
                                              index_col=0, \
                                              low_memory=True, error_bad_lines=True, \
                                              doublequote=False, memory_map=True, float_precision="high", \
                                              engine='c'
                                              )
            if not promfeatures:
                promfeatures = features
            else:
                promfeatures.merge(features.iloc[:, :-1], how='inner', left_index=True, right_index=True)
        labels = features.iloc[:, -1]

        """
        Making indexes equal for features and labels
        """
        promfeatures.merge(labels, how='inner', left_index=True, right_index=True)

        """
        Grabbing the source_sequence_id, features and labels
        """
        source_sequence_id = promfeatures.index
        labels = promfeatures.iloc[:, -1]
        promfeatures = promfeatures.iloc[:, :-1]

        assert promfeatures.shape[1] == 672

        return source_sequence_id, promfeatures, labels

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
        prop_as_feature = np.array(self.supInfo[prop_name]).astype(float)[seq_into_oligonucs_index]
        return prop_as_feature

    def properties_values_along_seqs(self,seq_list,prop_name_list):
        seqs_encoded_with_props = []
        for p in prop_name_list:
            seqs_encoded = []
            for s in seq_list:
                seqs_encoded.append(self.property_values_along_the_sequence(s,p))
            seqs_encoded_with_props.append(np.array(seqs_encoded))
        return np.array(seqs_encoded_with_props)

    def encode_fasta_into_pseknc_physicochemical_properties(self, input_path, prop_name,L=58, append_seq_id=True,
                                                            verbose=False):
        file_in = open(input_path,"r")
        sigma = ''
        feature_vector = []

        for i,line in enumerate(file_in):
            line_encoded = []
            if line == "":
                continue
            elif line[0] == '>':
                sigma = line[1:-1]
            else:
                if line[-1:] == "\n":
                    line58=line[:-1]
                else:
                    line58 = line

                if len(line58) != L:
                    print("Line {} has size {} instead of {}.".format(i,len(line58),L))
                    continue

                line_encoded = list(self.property_values_along_the_sequence(line58, prop_name))
                if verbose:
                    print("length da linha {}".format(len(line58)))
                    print("length da line_encoded {}".format(len(line_encoded)))

                if append_seq_id:
                    line_encoded.append(sigma)

            if line_encoded != []:
                feature_vector.append(line_encoded)

        file_in.close()

        return feature_vector
