import numpy as np
from git.beiko_lab.DNA_encoders.my_pseknc import Pseknc
import os
import pandas as pd


class Pc3mer(Pseknc):
    def __init__(self, classes=None, binary=True, folder_for_output=None,
                 file_name_prefix="", L=58):
        """

        :param classes:
        :param binary: True if we want the classes to be mapped into promoter or non-promoter only.
        False if we want to keep the promoter-type information, such as Sigma24, Sigma70, etc.
        Default: True
        :param folder_for_output:
        """
        # if binary is True and the parameter classes is not None,
        # the size of classes must be 2
        if binary and classes:
            assert len(classes) == 2

        super().__init__(binary=binary, L=L)

        self.name = 'pc3mer'
        self.folder_for_output = folder_for_output
        # if self.folder_for_output is not None:
        #     self.folder_for_output = os.path.join(self.folder_for_output, self.name)
        self.binary = binary
        self.file_name_prefix = file_name_prefix

        if classes:
            self.classes = np.array([self.synonyms(c) for c in classes])
            if binary:
                assert len(self.classes) == 2
        else:
            self.classes = None

    def property_values_along_seq(self, seq, prop_name):
        """
        Finds all the 3-mers that compose the sequence keeping the order they appeared in the sequence
        Then, finds the property's values for each one of these 3-mers
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

    def encode_fasta_with_prop(
        self, file_in, prop_name, L=58, append_seq_id=True, verbose=False
    ):
        """
        Receive the fasta file read with fid.readlines().
        For each pair of (">"+sequence_id+"\n"+sequence+"\n"), it converts sequence into pc3mer using
        only the physicochemical property named "prop_name" and search for a synonym for sequence_id among
        the possible classes.
        Append the encoded sequence and its class in the variable feature_vector and return feature_vector.
        :param file_in: Fasta file read using readlines()
        :param prop_name: Physicochemical property to use to encode the file.
        :param L: Sequence length.
        :param append_seq_id: Boolean. If True, it will append the class to the encoded sequence.
        :param verbose: Boolean. If True, it prints intermediate information.
        :return:
        """
        print(f"Encoding with property {prop_name} ...")

        sigma = ""
        feature_vector = []

        for i, line in enumerate(file_in):
            line_encoded = []
            if line == "":
                continue
            elif line.startswith(">"):
                sigma = line[1:-1]
            else:
                if line.endswith("\n"):
                    line58 = line[:-1]
                else:
                    line58 = line

                if len(line58) != L:
                    print(
                        "Line {} has size {} instead of {} and will not be encoded.".format(i, len(line58), L)
                    )
                    continue

                line_encoded = list(
                    self.property_values_along_seq(line58, prop_name)
                )
                if verbose:
                    print("line length {}".format(len(line58)))
                    print("encoded line length {}".format(len(line_encoded)))

                if append_seq_id:
                    line_encoded.append(self.synonyms(sigma))

            if line_encoded != []:
                feature_vector.append(line_encoded)

        return np.array(feature_vector)




    def encode_fasta_file(self, input_file_path, output_path=None, folder_for_output=None, L=58, props_names=None,
                          add_class_to_entries=True, verbose=False,
                          store_encode_by_indiv_prop=False):
        """
        Encodes each 3-mer in the sequence into physicochemical properties, store them in disk and return them.

        Input:
         :param input_file_path: (str) It is the path with the fasta file name and extension. The file is the one
         with DNA sequences that we want to encode into pc3mer
         :param output_path: (str) path and file name with extension. File where to save the fasta file encoded into
         one single pc3mer file. <Required> to obtain one single file for the fasta file encoded with all properties
         in props_names.
         :param folder_for_output: (str) It is the path where to create the Pc3mer folder with one file for each
         property used to encode the fasta file
         :param L: (int) Sequence length
         :param props_names: (list of str) List of properties to encode the sequence.
         Default: None. If None, all properties are used.
         :param add_class_to_entries: (bool) Add sequence's label to its encoded output.
         :param verbose: (bool) Print or not the intermediate results.
         :param store_encode_by_indiv_prop: (bool) Weather we want the fasta file to be encoded into one file per
         property or not. Set to <True> to obtain one file per property used for encoding.
        """

        print(f"Encoding {input_file_path} into pc3mer...")
        file_in = Pc3mer.read_fasta_file(input_file_path, L=L)

        # path to save outputs
        if folder_for_output is None:
            folder_for_output = self.folder_for_output

        if store_encode_by_indiv_prop:
            folder_for_output = os.path.join(folder_for_output, "Pc3mer")

        # if folder for output is still None, the encoded sequence will be returned
        # instead of saved in disk (see block of code at the end)
        if folder_for_output is not None:
            os.makedirs(folder_for_output, exist_ok=True)

        # encoding using the properties passed as parameter. Default: all properties
        if props_names is None:
            props_names = self.properties

        # create the df with the features of the first property
        feature_vector = pd.DataFrame(
            self.encode_fasta_with_prop(
                file_in, props_names[0], L=L, append_seq_id=add_class_to_entries, verbose=verbose
            )
        )
        # if the label was added to the end of feature_vector,
        # it is removed and stored to be used later
        if add_class_to_entries:
            labels = feature_vector.iloc[:, -1]
            feature_vector = feature_vector.iloc[:, :-1]

        # add the encoded sequences from the other properties to the df
        for prop_name in props_names[1:]:
            # encoding the whole file using property prop_name
            features_and_maybe_label = pd.DataFrame(
                self.encode_fasta_with_prop(
                    file_in, prop_name, L, append_seq_id=add_class_to_entries, verbose=verbose
                )
            )

            # storing encoded file per property
            if store_encode_by_indiv_prop:
                features_and_maybe_label.to_csv(
                    os.path.join(folder_for_output, self.file_name_prefix + prop_name + ".md"), header=False, index=False
                )

            # if the label is added to the end of the encoded sequence,
            # we need to remove before concatenating it
            if add_class_to_entries:
                labels_ = features_and_maybe_label.iloc[:, -1]
                feature_vector_ = features_and_maybe_label.iloc[:, :-1]
                feature_vector = pd.concat(
                    [
                        feature_vector,
                        feature_vector_
                    ], axis=1, ignore_index=True
                )
                # verify that current labels are exactly the same as the labels
                # from the first encoding
                assert labels.values.tolist() == labels_.values.tolist()
            else:
                feature_vector = pd.concat(
                    [
                        feature_vector,
                        features_and_maybe_label
                    ],
                    axis=1,
                    ignore_index=True,
                )

        # adding the labels to the final concatenation of the encoded seqs
        if add_class_to_entries:
            feature_vector = pd.concat(
                [feature_vector, labels], axis=1, ignore_index=True
            )

        if verbose:
            print(
                f"Done. {feature_vector.shape[0]} samples with {feature_vector.shape[1]-1} features."+\
                f"The statement 'label is at the end of each entry' is {add_class_to_entries}"
            )

        # save the df in disk
        if output_path is not None:
            feature_vector.to_csv(output_path, header=False, index=False)
            if verbose:
                print(f"Encoded sequence saved at {output_path}")
        elif folder_for_output is not None:
            file_out = os.path.join(folder_for_output, self.file_name_prefix+self.name+".md")
            feature_vector.to_csv(file_out, header=False, index=False)
            if verbose:
                print(f"Encoded sequence saved at {file_out}")

        return feature_vector

    @staticmethod
    def read_fasta_file(input_file_path, L=58):
        """
        Reads the file and verify if it is a valid fasta format file before
        returning it as a list of lines
        :param input_file_path:
        :return:
        """
        assert os.path.isfile(input_file_path), "The input fasta file {} does not exist.".format(input_file_path)

        with open(input_file_path, "r") as fid:
            file_in = fid.readlines()

        line0 = file_in[0]
        line1 = file_in[1][:-1]
        fid.close()
        if not(line0.startswith(">")) or len(line1) != L:
            print(
                f"The file at {input_file_path} is in the wrong format. Please, change it into fasta format with"+\
                f" sequences of {L}bp"
            )
            assert False, f"The file at {input_file_path} is in wrong format. Please, change it into fasta format "+\
                    f"with sequences with {L}bp"
        else:
            return file_in

    def convert_fasta_file_to_pc3mer_stats(self, input_file_path, folder_for_output=None):
        """
        Receive the fasta file read with fid.readlines().
        Encodes each 3-mer in the sequence into physicochemical properties and calculates its statistics
        Input:
         *fasta_file* is the path with the fasta file name and extension.
         *base_database* is the path where to create the folder with the stats files
        """
        print(f"Encoding {input_file_path} into pc3mer_stats...")
        file_in = Pc3mer.read_fasta_file(input_file_path)

        if folder_for_output is None:
            folder_for_output = self.folder_for_output

        folder_for_output = os.path.join(folder_for_output, "Pc3mer_stats")
        os.makedirs(folder_for_output, exist_ok=True)

        for prop in self.properties:
            pse = pd.DataFrame(
                self.encode_fasta_with_prop(
                    file_in, prop
                )
            )
            pse_stats = pd.DataFrame([])
            sigma = pse.iloc[:, 56].values
            temp = pse.iloc[:, :56].astype('float').values

            pse_stats["min"] = np.min(temp, axis=1)
            pse_stats["max"] = np.max(temp, axis=1)
            pse_stats["mean"] = np.mean(temp, axis=1)
            pse_stats["std"] = np.std(temp, axis=1)
            pse_stats["median"] = np.median(temp, axis=1)
            pse_stats["variance"] = np.var(temp, axis=1)
            pse_stats["Sigma"] = sigma
            del temp
            del pse

            pse_stats = pd.DataFrame(data=pse_stats)
            pse_stats.to_csv(
                os.path.join(folder_for_output, prop + ".md"), header=False, index=False
            )
            del pse_stats