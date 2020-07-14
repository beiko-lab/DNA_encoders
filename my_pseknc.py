import numpy as np

class pseknc():

    def __init__(self,L=58):

        self.oligonucs = ["GGG","GGA","GGC","GGT","GAG","GAA","GAC","GAT","GCG","GCA","GCC","GCT","GTG","GTA","GTC","GTT","AGG","AGA",\
           "AGC","AGT","AAG","AAA","AAC","AAT","ACG","ACA","ACC","ACT","ATG","ATA","ATC","ATT","CGG","CGA","CGC","CGT",
           "CAG","CAA","CAC","CAT","CCG","CCA","CCC","CCT","CTG","CTA","CTC","CTT","TGG","TGA","TGC","TGT","TAG","TAA",
           "TAC","TAT","TCG","TCA","TCC","TCT","TTG","TTA","TTC","TTT"]
        self.supInfo = {"Bendability-DNAse":[5.7,6.2,8.2,5.2,6.6,5.1,5.6,3.6,4.3,7.5,8.2,6.3,6.8,6.4,5.6,1.6,4.7,6.5,6.3,2,4.2,0.1,\
                                1.6,0,5.2,5.8,5.2,2,8.7,9.7,3.6,0,3,5.8,4.3,5.2,9.6,6.2,6.8,8.7,3,0.7,5.7,4.7,9.6,7.8,\
                                6.6,4.2,0.7,10,7.5,5.8,7.8,7.3,6.4,9.7,5.8,10,6.2,6.5,6.2,7.3,5.1,0.1],\
           "Bendability-consensus": [5.85,5,9.1,5.3,6,4.05,5.5,4.45,5.9,6.75,9.1,6.9,6.65,5.05,5.5,2.65,5.05,4.9,6.9,
                                     3.9,4.7,0.05,2.65,0.35,5.3,5.5,5.3,3.9,7.7,6.25,4.45,0.35,3.85,7.05,5.9,5.3,6.9,
                                     4.75,6.65,7.7,3.85,3.05,5.85,5.05,6.9,5,6,4.7,3.05,7.7,6.75,5.5,5,4.65,5.05,6.25,
                                     7.05,7.7,5,4.9,4.75,4.65,4.05,0.05],\
           "Trinucleotide GC Content": [3,2,3,2,2,1,2,1,3,2,3,2,2,1,2,1,2,1,2,1,1,0,1,0,2,1,2,1,1,0,1,0,3,2,3,2,2,1,2,1,\
                                        3,2,3,2,2,1,2,1,2,1,2,1,1,0,1,0,2,1,2,1,1,0,1,0],\
           "Nucleosome positioning": [13,-5,45,8,8,-12,8,7,25,13,45,25,17,-6,8,-6,8,-9,25,11,6,-36,-6,-30,8,6,8,11,18,\
                                      -13,7,-30,2,31,25,8,-2,-9,17,18,2,8,13,8,-2,-18,8,6,8,8,13,6,-18,-20,-6,-13,31,8,\
                                      -5,-9,-9,-20,-12,-36],\
           "Consensus_roll": [5.827,4.9907,9.0823,5.31605,5.9806,4.0633,5.51645,4.44325,5.89135,6.75525,9.0823,6.8829,\
                              6.62555,5.0673,5.51645,2.64115,5.0523,4.8884,6.8829,3.9232,4.69915,0.0633,2.64115,0.35,\
                              5.3055,5.4903,5.31605,3.9232,7.7171,6.2734,4.44325,0.35,3.869,7.07195,5.89135,5.3055,\
                              6.8996,4.7618,6.62555,7.7171,3.869,3.05865,5.827,5.0523,6.8996,5.00295,5.9806,4.69915,\
                              3.05865,7.7,6.75525,5.4903,5.00295,4.6709,5.0673,6.2734,7.07195,7.7,4.9907,4.8884,4.7618,\
                              4.6709,4.0633,0.0633],\
           "Consensus_Rigid": [5.827,4.9907,9.0823,5.31605,5.9806,4.0633,5.51645,4.44325,5.89135,6.75525,9.0823,6.8829,\
                               6.62555,5.0673,5.51645,2.64115,5.0523,4.8884,6.8829,3.9232,4.69915,0.0633,2.64115,0.35,\
                               5.3055,5.4903,5.31605,3.9232,7.7171,6.2734,4.44325,0.35,3.869,7.07195,5.89135,5.3055,\
                               6.8996,4.7618,6.62555,7.7171,3.869,3.05865,5.827,5.0523,6.8996,5.00295,5.9806,4.69915,\
                               3.05865,7.7,6.75525,5.4903,5.00295,4.6709,5.0673,6.2734,7.07195,7.7,4.9907,4.8884,4.7618,\
                               4.6709,4.0633,0.0633],\
           "Dnase I": [3.311,3.819,1.387,3.619,3.221,4.385,3.498,4.153,3.275,2.754,1.387,2.683,2.832,3.77,3.498,5.26,\
                       3.782,3.879,2.683,4.471,3.995,6.882,5.26,6.698,3.625,3.516,3.619,4.471,2.185,3.047,4.153,6.698,\
                       4.502,2.57,3.275,3.625,2.671,3.958,2.832,2.185,4.502,5,3.311,3.782,2.671,3.813,3.221,3.995,5,10,\
                       2.754,3.516,3.813,4.013,3.77,3.047,2.57,2.197,3.819,3.879,3.958,4.013,4.385,0.1],\
           "Dnase I-Rigid": [3.868,3.581,2.448,4.156,3.353,4.214,3.925,5.087,4.678,2.842,2.448,3.524,3.239,3.467,3.925,\
                             6.272,4.445,3.41,3.524,6.033,4.736,7.176,6.272,7.237,4.156,3.81,4.156,6.033,2.169,1.613,\
                             5.087,7.237,5.44,3.81,4.678,4.156,1.668,3.581,3.239,2.169,5.44,6.813,3.868,4.445,1.668,\
                             2.673,3.353,4.736,6.813,1.447,2.842,3.81,2.673,2.955,3.467,1.613,3.81,1.447,3.581,3.41,\
                             3.581,2.955,4.214,7.176],\
           "MW-Daltons": [622.4,622.4,622.4,622.4,621.4,621.4,621.4,621.4,622.4,622.4,622.4,622.4,621.4,621.4,621.4,\
                          621.4,622.4,622.4,622.4,622.4,621.4,621.4,621.4,621.4,622.4,622.4,622.4,622.4,621.4,621.4,\
                          621.4,621.4,622.4,622.4,622.4,622.4,621.4,621.4,621.4,621.4,622.4,622.4,622.4,622.4,621.4,\
                          621.4,621.4,621.4,622.4,622.4,622.4,622.4,621.4,621.4,621.4,621.4,622.4,622.4,622.4,622.4,\
                          621.4,621.4,621.4,621.4],\
           "MW-kg": [103.3887,103.3887,103.3887,103.3887,103.22259,103.22259,103.22259,103.22259,103.3887,103.3887,\
                     103.3887,103.3887,103.22259,103.22259,103.22259,103.22259,103.3887,103.3887,103.3887,103.3887,\
                     103.22259,103.22259,103.22259,103.22259,103.3887,103.3887,103.3887,103.3887,103.22259,103.22259,\
                     103.22259,103.22259,103.3887,103.3887,103.3887,103.3887,103.22259,103.22259,103.22259,103.22259,\
                     103.3887,103.3887,103.3887,103.3887,103.22259,103.22259,103.22259,103.22259,103.3887,103.3887,\
                     103.3887,103.3887,103.22259,103.22259,103.22259,103.22259,103.3887,103.3887,103.3887,103.3887,\
                     103.22259,103.22259,103.22259,103.22259],\
           "Nucleosome": [6,3.8,10,5.4,5.4,3,5.4,5.3,7.5,6,10,7.5,6.5,3.7,5.4,3.7,5.4,3.3,7.5,5.8,5.2,0,3.7,0.7,5.4,\
                          5.2,5.4,5.8,6.7,2.8,5.3,0.7,4.7,8.3,7.5,5.4,4.2,3.3,6.5,6.7,4.7,5.4,6,5.4,4.2,2.2,5.4,5.2,\
                          5.4,5.4,6,5.2,2.2,2,3.7,2.8,8.3,5.4,3.8,3.3,3.3,2,3,0],\
           "Nucleosome-Rigid": [3.536,4.799,1.309,3.878,3.878,5.264,3.878,3.935,2.691,3.536,1.309,2.691,3.253,4.857,\
                                3.878,4.857,3.878,5.089,2.691,3.65,3.992,7.045,4.857,6.624,3.878,3.992,3.878,3.65,3.14,\
                                5.381,3.935,6.624,4.279,2.245,2.691,3.878,4.567,5.089,3.253,3.14,4.279,3.878,3.536,\
                                3.878,4.567,5.734,3.878,3.992,3.878,3.878,3.536,3.992,5.734,5.852,4.857,5.381,2.245,\
                                3.878,4.799,5.089,5.089,5.852,5.264,7.045]
           }
        self.properties = ['Bendability-DNAse', 'Bendability-consensus', 'Trinucleotide GC Content', 'Nucleosome positioning',
         'Consensus_roll', 'Consensus_Rigid', 'Dnase I', 'Dnase I-Rigid', 'MW-Daltons', 'MW-kg', 'Nucleosome',
         'Nucleosome-Rigid']
        self.lamb = 2
        self.L = L
        self.p_v()
        self.w = 0.3

    def p_v(self):
        """
        pv refers to P_v, which is strutural property v
        P_v(nucleotide1,nucleotide2,nucleotide3) =
            (P_v(nucleotide1,nucleotide2,nucleotide3) - mean(P_v for all 3-tuple nucleotides))
            /
            std(P_v for all 3-tuple nucleotides)

        :return: supInfo normalized
        """
        for props in self.properties:
            originalValues = self.supInfo[props]
            _mean = np.mean(originalValues)
            _std = np.std(originalValues)
            self.supInfo[props] = (originalValues - _mean ) / _std

    def encode_into_kmers(self,seq):
        """Count kmer occurrences in a given seq.
        Parameters
        ----------
        seq : string - A single DNA sequence.
        Returns
        -------
        seq_into_oligonucs_index : list of size len(seq)-2 with the oligonucs' indexes of each 3-mer

        DataFrame, columns are oligonucleotides and line 0 has the # of times each oligonucleotide appears in seq
            A dictionary of counts keyed by their individual kmers (strings
            of length k).
        """
        seq_into_oligonucs_index = []
        kmerCounts = np.zeros([len(self.oligonucs)])

        # loop over the sequence grabbing all 3-mers
        for l in range(len(seq)-2):
            kmer = seq[l:l+3]
            try:
                # searches oligonucs for the specific kmer and return its index in that array
                index = self.oligonucs.index(kmer)
            except:
                print("The kmer {} is not in the possible oligonucs {}".format(kmer,self.oligonucs))
            # append the kmer index number to create a list of kmers indexes (these indexes are in the oligonucs)
            seq_into_oligonucs_index.append(index)
            # add 1 to the count of what and how many kmers were found
            kmerCounts[index] += 1

        # find the normalized frequence for the kmers in the sequence
        normal_freq_of_kmers = kmerCounts/np.sum(kmerCounts)

        return seq_into_oligonucs_index, normal_freq_of_kmers

    def big_theta(self,index_kmer1,index_kmer2):
        ni = len(self.properties) #ni = number of properties being extracted
        props_sum_for_kmers_at_i_and_ij = 0
        # sum of (p_v(kmer at pos i) - p_v(kmer at pos i+j))^2 for v in properties
        for v in self.properties:
            props_sum_for_kmers_at_i_and_ij += np.power(self.supInfo[v][index_kmer1] - self.supInfo[v][index_kmer2],2)

        # devide the sum for ni = number of properties being extracted
        return props_sum_for_kmers_at_i_and_ij/ni

    def small_theta(self,seq_into_oligonucs_index):
        #[0, 1, 5, 22, 24, 35, 12, 51, 13]

        theta = np.zeros([self.lamb])

        for j in range(1,self.lamb+1):
            # starts the usm with zero and adding the big_theta results
            big_theta_sum=0
            for i in range(1,(self.L-j-1)):
                #subtracting 1 to access the index of array seq_into_oligonucs_index
                #print('\n')
                #print(len(seq_into_oligonucs_index))
                #print(i-1)
                #print(i+j-1)
                big_theta_sum += self.big_theta(seq_into_oligonucs_index[i-1],seq_into_oligonucs_index[i+j-1])

            # j-1 is used to access since position 1 (equ. index=0) of the array
            theta[j-1] = big_theta_sum/(self.L-j-1)

        return theta

    def encode_into_pseknc(self,seq):
        seq_into_oligonucs_index, normal_freq_of_kmers = self.encode_into_kmers(seq)
        # w_theta is w * small_theta
        #print("len of seq {}".format(len(seq)))
        #print("len of seq_into_oligonucs_index {}".format(len(seq_into_oligonucs_index)))
        w_theta = self.w * self.small_theta(seq_into_oligonucs_index)
        feature_vector = (np.concatenate((normal_freq_of_kmers,w_theta))/(np.sum(normal_freq_of_kmers)+np.sum(w_theta)))
        return feature_vector

    def encode_fasta_into_pseknc(self, input_path, output_path, save = True):
        file_in = open(input_path,"r")
        if save:
            file_out = open(output_path, "w")
        else:
            encoded = []
        sigma = ''

        for i,line in enumerate(file_in):
            if (line == ""):
                continue
            elif(line[0] == '>'):
                sigma = line[1:-1]
            else:
                if(line[-1:] == "\n"):
                    line58=line[:-1]
                else:
                    line58 = line

                if(len(line58) != self.L):
                    print("Line {} has size {} instead of {}.".format(i,len(line58),self.L))
                    continue
                line_enconded = self.encode_into_pseknc(line58)
                if save:
                    file_out.write(str(list(line_enconded))[1:-1]+","+sigma+"\n")
                else:
                    encoded.append(line_enconded)
        file_in.close()
        if save:
            file_out.close()
            return output_path
        else:
            return encoded


