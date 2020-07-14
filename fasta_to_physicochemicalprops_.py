import pandas as pd
import os.path as osp
import numpy as np
import os


"""
This encoder reads a fasta file and encodes the sequences' 3-mers into 12 files, one for each physicochemical property.

The fasta file needs to have the format:
        >1
        GGTTTATTGCCTTGCAGCTGGCGAGAGACGGTATTGCTCATGCACAAGCCTTGTTCAG
        >2
        TGCCCTGACTTCACCCCGCTGTGTCTGCTTTTCCCGACTATTCTTAATGAGCTTCGAT
        >3
        AATGTGGATAGATATGAATTATTTTTCTCCTTAAGGATCATCCGTTATTTGGGTCGTT
        >4
        CAGTTATTTACCTTACTTTACGCGCGCGTAACTCTGGCAACATCACTACAGGATAGCG
        >5
        AAAAAGTTATACGCGGTGGAAACATTGCCCGGATAGTCTATAGTCACTAAGCATTAAA
"""

from physicochemicalfor3mer import PhysicochemicalFor3mer


class FastaToPhysicochemical3mers():

    def __init__(self, folder_for_output=None):
        # creates the encoder
        self.folder_for_output=folder_for_output
        self.encoder = PhysicochemicalFor3mer(folder_for_output=folder_for_output)
        self.props = self.encoder.properties

    def convert_to_physicochemical_properties(self, input_file_path, folder_for_output=None):
        """
        Encodes each 3-mer in the sequence into physicochemical properties
        Input:
         *fasta_file* is the path with the fasta file name and extension.
         *base_database* is the path where to create the folder with the properties files
        """
        self.secure_input_fasta(input_file_path)

        if folder_for_output is None:
            folder_for_output = self.folder_for_output

        if not osp.isdir(folder_for_output + "Physicochemical_props/"):
            os.mkdir(folder_for_output + "Physicochemical_props/")

        for prop in self.props:
            pse = pd.DataFrame(self.encoder.encode_fasta_into_pseknc_physicochemical_properties(input_file_path, prop))
            pse.to_csv(folder_for_output + "Physicochemical_props/" + prop + ".md",
                       header=False, index=False)
            del pse

    def secure_input_fasta(self,input_file_path):
        if not osp.isfile(input_file_path):
            print("The input fasta file {} does not exist.".format(input_file_path))
            assert False
        else:
            fid = open(input_file_path,"r")
            line0 = fid.readline()
            line1 = fid.readline()[:-1]
            fid.close()
            if line0[0] != ">" or len(line1) != 58:
                print("The file at {} is in wrong format. Please, change it into fasta format with sequences with 58bp".
                      format(input_file_path))
                assert False

    def convert_to_physicochemical_stats(self, input_file_path, folder_for_output=None):
        """
        Encodes each 3-mer in the sequence into physicochemical properties and calculates its statistics
        Input:
         *fasta_file* is the path with the fasta file name and extension.
         *base_database* is the path where to create the folder with the stats files
        """
        self.secure_input_fasta(input_file_path)

        if folder_for_output is None:
            folder_for_output = self.folder_for_output

        if not osp.isdir(folder_for_output + "Physicochemical_stats/"):
            os.mkdir(folder_for_output + "Physicochemical_stats/")

        for prop in self.props:
            pse = pd.DataFrame(self.encoder.encode_fasta_into_pseknc_physicochemical_properties(input_file_path,prop))
            pse_stats = pd.DataFrame([])
            sigma = pse.iloc[:,56].values
            temp = pse.iloc[:,:56].values

            pse_stats["min"] = np.min(temp,axis=1)
            pse_stats["max"] = np.max(temp,axis=1)
            pse_stats["mean"] = np.mean(temp,axis=1)
            pse_stats["std"] = np.std(temp,axis=1)
            pse_stats["median"] = np.median(temp,axis=1)
            pse_stats["variance"] = np.var(temp,axis=1)
            pse_stats["Sigma"] = sigma
            del temp
            del pse

            pse_stats = pd.DataFrame(data=pse_stats)
            pse_stats.to_csv(folder_for_output + "Physicochemical_stats/" + prop + ".md",
                             header = False, index = False)
            del pse_stats







