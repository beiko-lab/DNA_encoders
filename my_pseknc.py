import numpy as np
from git.beiko_lab.DNA_encoders.base_for_encoders import BaseForEncoder

class Pseknc(BaseForEncoder):

    def __init__(self, L=58, binary=True):
        super().__init__(binary=binary)
        self.oligonucs = ["GGG", "GGA", "GGC", "GGT", "GAG", "GAA", "GAC", "GAT", "GCG", "GCA", "GCC", "GCT", "GTG",
                          "GTA", "GTC", "GTT", "AGG", "AGA", "AGC", "AGT", "AAG", "AAA", "AAC", "AAT", "ACG", "ACA",
                          "ACC", "ACT", "ATG", "ATA", "ATC", "ATT", "CGG", "CGA", "CGC", "CGT", "CAG", "CAA", "CAC",
                          "CAT", "CCG", "CCA", "CCC", "CCT", "CTG", "CTA", "CTC", "CTT", "TGG", "TGA", "TGC", "TGT",
                          "TAG", "TAA", "TAC", "TAT", "TCG", "TCA", "TCC", "TCT", "TTG", "TTA", "TTC", "TTT"]
        self.original_supInfo = {
           "Bendability-DNAse": [5.7, 6.2, 8.2, 5.2, 6.6, 5.1, 5.6, 3.6, 4.3, 7.5, 8.2, 6.3, 6.8, 6.4, 5.6, 1.6,
                                 4.7, 6.5, 6.3, 2, 4.2, 0.1,  1.6, 0, 5.2, 5.8, 5.2, 2, 8.7, 9.7, 3.6, 0, 3, 5.8, 4.3,
                                 5.2, 9.6, 6.2, 6.8, 8.7, 3, 0.7, 5.7, 4.7, 9.6, 7.8, 6.6, 4.2, 0.7, 10, 7.5, 5.8,
                                 7.8, 7.3, 6.4, 9.7, 5.8, 10, 6.2, 6.5, 6.2, 7.3, 5.1, 0.1],
           "Bendability-consensus": [5.85, 5, 9.1, 5.3, 6, 4.05, 5.5, 4.45, 5.9, 6.75, 9.1, 6.9, 6.65, 5.05, 5.5, 2.65,
                                     5.05, 4.9, 6.9, 3.9, 4.7, 0.05, 2.65, 0.35, 5.3, 5.5, 5.3, 3.9, 7.7, 6.25, 4.45,
                                     0.35, 3.85, 7.05, 5.9, 5.3, 6.9, 4.75, 6.65, 7.7, 3.85, 3.05, 5.85, 5.05, 6.9, 5,
                                     6, 4.7, 3.05, 7.7, 6.75, 5.5, 5, 4.65, 5.05, 6.25, 7.05, 7.7, 5, 4.9, 4.75, 4.65,
                                     4.05, 0.05],
           "Trinucleotide GC Content": [3, 2, 3, 2, 2, 1, 2, 1, 3, 2, 3, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 0, 1, 0, 2, 1, 2,
                                        1, 1, 0, 1, 0, 3, 2, 3, 2, 2, 1, 2, 1, 3, 2, 3, 2, 2, 1, 2, 1, 2, 1, 2, 1, 1, 0,
                                        1, 0, 2, 1, 2, 1, 1, 0, 1, 0],
           "Nucleosome positioning": [13, -5, 45, 8, 8, -12, 8, 7, 25, 13, 45, 25, 17, -6, 8, -6, 8, -9, 25, 11, 6, -36,
                                      -6, -30, 8, 6, 8, 11, 18, -13, 7, -30, 2, 31, 25, 8, -2, -9, 17, 18, 2, 8, 13, 8,
                                      -2, -18, 8, 6, 8, 8, 13, 6, -18, -20, -6, -13, 31, 8, -5, -9, -9, -20, -12, -36],
           "Consensus_roll": [5.827, 4.9907, 9.0823, 5.31605, 5.9806, 4.0633, 5.51645, 4.44325, 5.89135, 6.75525,
                              9.0823, 6.8829, 6.62555, 5.0673, 5.51645, 2.64115, 5.0523, 4.8884, 6.8829, 3.9232,
                              4.69915, 0.0633, 2.64115, 0.35, 5.3055, 5.4903, 5.31605, 3.9232, 7.7171, 6.2734, 4.44325,
                              0.35, 3.869, 7.07195, 5.89135, 5.3055, 6.8996, 4.7618, 6.62555, 7.7171, 3.869, 3.05865,
                              5.827, 5.0523, 6.8996, 5.00295, 5.9806, 4.69915, 3.05865, 7.7, 6.75525, 5.4903, 5.00295,
                              4.6709, 5.0673, 6.2734, 7.07195, 7.7, 4.9907, 4.8884, 4.7618, 4.6709, 4.0633, 0.0633],
           "Consensus_Rigid": [5.827, 4.9907, 9.0823, 5.31605, 5.9806, 4.0633, 5.51645, 4.44325, 5.89135, 6.75525,
                               9.0823, 6.8829, 6.62555, 5.0673, 5.51645, 2.64115, 5.0523, 4.8884, 6.8829, 3.9232,
                               4.69915, 0.0633, 2.64115, 0.35, 5.3055, 5.4903, 5.31605, 3.9232, 7.7171, 6.2734, 4.44325,
                               0.35, 3.869, 7.07195, 5.89135, 5.3055, 6.8996, 4.7618, 6.62555, 7.7171, 3.869, 3.05865,
                               5.827, 5.0523, 6.8996, 5.00295, 5.9806, 4.69915, 3.05865, 7.7, 6.75525, 5.4903, 5.00295,
                               4.6709, 5.0673, 6.2734, 7.07195, 7.7, 4.9907, 4.8884, 4.7618, 4.6709, 4.0633, 0.0633],
           "Dnase I": [3.311, 3.819, 1.387, 3.619, 3.221, 4.385, 3.498, 4.153, 3.275, 2.754, 1.387, 2.683, 2.832, 3.77,
                       3.498, 5.26, 3.782, 3.879, 2.683, 4.471, 3.995, 6.882, 5.26, 6.698, 3.625, 3.516, 3.619, 4.471,
                       2.185, 3.047, 4.153, 6.698, 4.502, 2.57, 3.275, 3.625, 2.671, 3.958, 2.832, 2.185, 4.502, 5,
                       3.311, 3.782, 2.671, 3.813, 3.221, 3.995, 5, 10, 2.754, 3.516, 3.813, 4.013, 3.77, 3.047, 2.57,
                       2.197, 3.819, 3.879, 3.958, 4.013, 4.385, 0.1],
           "Dnase I-Rigid": [3.868, 3.581, 2.448, 4.156, 3.353, 4.214, 3.925, 5.087, 4.678, 2.842, 2.448, 3.524, 3.239,
                             3.467, 3.925, 6.272, 4.445, 3.41, 3.524, 6.033, 4.736, 7.176, 6.272, 7.237, 4.156, 3.81,
                             4.156, 6.033, 2.169, 1.613, 5.087, 7.237, 5.44, 3.81, 4.678, 4.156, 1.668, 3.581, 3.239,
                             2.169, 5.44, 6.813, 3.868, 4.445, 1.668, 2.673, 3.353, 4.736, 6.813, 1.447, 2.842, 3.81,
                             2.673, 2.955, 3.467, 1.613, 3.81, 1.447, 3.581, 3.41, 3.581, 2.955, 4.214, 7.176],
           "MW-Daltons": [622.4, 622.4, 622.4, 622.4, 621.4, 621.4, 621.4, 621.4, 622.4, 622.4, 622.4, 622.4, 621.4,
                          621.4, 621.4, 621.4, 622.4, 622.4, 622.4, 622.4, 621.4, 621.4, 621.4, 621.4, 622.4, 622.4,
                          622.4, 622.4, 621.4, 621.4, 621.4, 621.4, 622.4, 622.4, 622.4, 622.4, 621.4, 621.4, 621.4,
                          621.4, 622.4, 622.4, 622.4, 622.4, 621.4, 621.4, 621.4, 621.4, 622.4, 622.4, 622.4, 622.4,
                          621.4, 621.4, 621.4, 621.4, 622.4, 622.4, 622.4, 622.4, 621.4, 621.4, 621.4, 621.4],
           "MW-kg": [103.3887, 103.3887, 103.3887, 103.3887, 103.22259, 103.22259, 103.22259, 103.22259, 103.3887,
                     103.3887, 103.3887, 103.3887, 103.22259, 103.22259, 103.22259, 103.22259, 103.3887, 103.3887,
                     103.3887, 103.3887, 103.22259, 103.22259, 103.22259, 103.22259, 103.3887, 103.3887, 103.3887,
                     103.3887, 103.22259, 103.22259, 103.22259, 103.22259, 103.3887, 103.3887, 103.3887, 103.3887,
                     103.22259, 103.22259, 103.22259, 103.22259, 103.3887, 103.3887, 103.3887, 103.3887, 103.22259,
                     103.22259, 103.22259, 103.22259, 103.3887, 103.3887, 103.3887, 103.3887, 103.22259, 103.22259,
                     103.22259, 103.22259, 103.3887, 103.3887, 103.3887, 103.3887, 103.22259, 103.22259, 103.22259,
                     103.22259],
           "Nucleosome": [6, 3.8, 10, 5.4, 5.4, 3, 5.4, 5.3, 7.5, 6, 10, 7.5, 6.5, 3.7, 5.4, 3.7, 5.4, 3.3, 7.5, 5.8,
                          5.2, 0, 3.7, 0.7, 5.4, 5.2, 5.4, 5.8, 6.7, 2.8, 5.3, 0.7, 4.7, 8.3, 7.5, 5.4, 4.2, 3.3, 6.5,
                          6.7, 4.7, 5.4, 6, 5.4, 4.2, 2.2, 5.4, 5.2, 5.4, 5.4, 6, 5.2, 2.2, 2, 3.7, 2.8, 8.3, 5.4, 3.8,
                          3.3, 3.3, 2, 3, 0],
           "Nucleosome-Rigid": [3.536, 4.799, 1.309, 3.878, 3.878, 5.264, 3.878, 3.935, 2.691, 3.536, 1.309, 2.691,
                                3.253, 4.857, 3.878, 4.857, 3.878, 5.089, 2.691, 3.65, 3.992, 7.045, 4.857, 6.624,
                                3.878, 3.992, 3.878, 3.65, 3.14, 5.381, 3.935, 6.624, 4.279, 2.245, 2.691, 3.878, 4.567,
                                5.089, 3.253, 3.14, 4.279, 3.878, 3.536, 3.878, 4.567, 5.734, 3.878, 3.992, 3.878,
                                3.878, 3.536, 3.992, 5.734, 5.852, 4.857, 5.381, 2.245, 3.878, 4.799, 5.089, 5.089,
                                5.852, 5.264, 7.045]
           }
        self.supInfo = {
         'Bendability-DNAse': [0.07230364, 0.26511335, 1.0363522, -0.12050607, 0.41936112,
                                     -0.15906801, 0.0337417, -0.73749715, -0.46756355, 0.7664186,
                                     1.0363522, 0.3036753, 0.49648501, 0.34223724, 0.0337417,
                                     -1.50873599, -0.31331578, 0.38079918, 0.3036753, -1.35448822,
                                     -0.50612549, -2.08716513, -1.50873599, -2.12572707, -0.12050607,
                                     0.11086558, -0.12050607, -1.35448822, 1.22916191, 1.61478134,
                                     -0.73749715, -2.12572707, -0.9688688, 0.11086558, -0.46756355,
                                     -0.12050607, 1.57621939, 0.26511335, 0.49648501, 1.22916191,
                                     -0.9688688, -1.85579348, 0.07230364, -0.31331578, 1.57621939,
                                     0.88210443, 0.41936112, -0.50612549, -1.85579348, 1.73046716,
                                     0.7664186, 0.11086558, 0.88210443, 0.68929472, 0.34223724,
                                     1.61478134, 0.11086558, 1.73046716, 0.26511335, 0.38079918,
                                     0.26511335, 0.68929472, -0.15906801, -2.08716513],
         'Bendability-consensus': [0.3577835, -0.09696936, 2.09654444, 0.06353165, 0.43803401,
                                         -0.60522256, 0.17053232, -0.39122121, 0.38453367, 0.83928653,
                                         2.09654444, 0.91953704, 0.7857862, -0.07021919, 0.17053232,
                                         -1.35422727, -0.07021919, -0.1504697, 0.91953704, -0.68547306,
                                         -0.25747037, -2.74523603, -1.35422727, -2.58473502, 0.06353165,
                                         0.17053232, 0.06353165, -0.68547306, 1.34753973, 0.57178485,
                                         -0.39122121, -2.58473502, -0.71222323, 0.99978754, 0.38453367,
                                         0.06353165, 0.91953704, -0.2307202, 0.7857862, 1.34753973,
                                         -0.71222323, -1.14022593, 0.3577835, -0.07021919, 0.91953704,
                                         -0.09696936, 0.43803401, -0.25747037, -1.14022593, 1.34753973,
                                         0.83928653, 0.17053232, -0.09696936, -0.28422054, -0.07021919,
                                         0.57178485, 0.99978754, 1.34753973, -0.09696936, -0.1504697,
                                         -0.2307202, -0.28422054, -0.60522256, -2.74523603],
         'Trinucleotide GC Content': [1.73205081, 0.57735027, 1.73205081, 0.57735027, 0.57735027,
                                            -0.57735027, 0.57735027, -0.57735027, 1.73205081, 0.57735027,
                                            1.73205081, 0.57735027, 0.57735027, -0.57735027, 0.57735027,
                                            -0.57735027, 0.57735027, -0.57735027, 0.57735027, -0.57735027,
                                            -0.57735027, -1.73205081, -0.57735027, -1.73205081, 0.57735027,
                                            -0.57735027, 0.57735027, -0.57735027, -0.57735027, -1.73205081,
                                            -0.57735027, -1.73205081, 1.73205081, 0.57735027, 1.73205081,
                                            0.57735027, 0.57735027, -0.57735027, 0.57735027, -0.57735027,
                                            1.73205081, 0.57735027, 1.73205081, 0.57735027, 0.57735027,
                                            -0.57735027, 0.57735027, -0.57735027, 0.57735027, -0.57735027,
                                            0.57735027, -0.57735027, -0.57735027, -1.73205081, -0.57735027,
                                            -1.73205081, 0.57735027, -0.57735027, 0.57735027, -0.57735027,
                                            -0.57735027, -1.73205081, -0.57735027, -1.73205081],
         'Nucleosome positioning': [0.57187906, -0.50109273, 2.47938447, 0.27383134, 0.27383134,
                                          -0.91835954, 0.27383134, 0.2142218, 1.28719359, 0.57187906,
                                          2.47938447, 1.28719359, 0.81031724, -0.56070227, 0.27383134,
                                          -0.56070227, 0.27383134, -0.7395309, 1.28719359, 0.45265997,
                                          0.15461225, -2.34898859, -0.56070227, -1.99133133, 0.27383134,
                                          0.15461225, 0.27383134, 0.45265997, 0.86992678, -0.97796908,
                                          0.2142218, -1.99133133, -0.08382592, 1.64485085, 1.28719359,
                                          0.27383134, -0.3222641, -0.7395309, 0.81031724, 0.86992678,
                                          -0.08382592, 0.27383134, 0.57187906, 0.27383134, -0.3222641,
                                          -1.2760168, 0.27383134, 0.15461225, 0.27383134, 0.27383134,
                                          0.57187906, 0.15461225, -1.2760168, -1.39523589, -0.56070227,
                                          -0.97796908, 1.64485085, 0.27383134, -0.50109273, -0.7395309,
                                          -0.7395309, -1.39523589, -0.91835954, -2.34898859],
         'Consensus_roll': [0.34471906, -0.10347006, 2.08929645, 0.07089121, 0.42703623,
                                  -0.60048142, 0.17828939, -0.39685897, 0.37920546, 0.84218595,
                                  2.08929645, 0.91059602, 0.77267725, -0.06241866, 0.17828939,
                                  -1.36263874, -0.07045745, -0.15829459, 0.91059602, -0.67556368,
                                  -0.25971727, -2.74415775, -1.36263874, -2.59050975, 0.06523726,
                                  0.16427511, 0.07089121, -0.67556368, 1.35765972, 0.58395334,
                                  -0.39685897, -2.59050975, -0.7046105, 1.01191152, 0.37920546,
                                  0.06523726, 0.91954587, -0.22614194, 0.77267725, 1.35765972,
                                  -0.7046105, -1.13889253, 0.34471906, -0.07045745, 0.91954587,
                                  -0.09690506, 0.42703623, -0.25971727, -1.13889253, 1.3484955,
                                  0.84218595, 0.16427511, -0.09690506, -0.27485699, -0.06241866,
                                  0.58395334, 1.01191152, 1.3484955, -0.10347006, -0.15829459,
                                  -0.22614194, -0.27485699, -0.60048142, -2.74415775],
         'Consensus_Rigid': [0.34471906, -0.10347006, 2.08929645, 0.07089121, 0.42703623,
                                   -0.60048142, 0.17828939, -0.39685897, 0.37920546, 0.84218595,
                                   2.08929645, 0.91059602, 0.77267725, -0.06241866, 0.17828939,
                                   -1.36263874, -0.07045745, -0.15829459, 0.91059602, -0.67556368,
                                   -0.25971727, -2.74415775, -1.36263874, -2.59050975, 0.06523726,
                                   0.16427511, 0.07089121, -0.67556368, 1.35765972, 0.58395334,
                                   -0.39685897, -2.59050975, -0.7046105, 1.01191152, 0.37920546,
                                   0.06523726, 0.91954587, -0.22614194, 0.77267725, 1.35765972,
                                   -0.7046105, -1.13889253, 0.34471906, -0.07045745, 0.91954587,
                                   -0.09690506, 0.42703623, -0.25971727, -1.13889253, 1.3484955,
                                   0.84218595, 0.16427511, -0.09690506, -0.27485699, -0.06241866,
                                   0.58395334, 1.01191152, 1.3484955, -0.10347006, -0.15829459,
                                   -0.22614194, -0.27485699, -0.60048142, -2.74415775],
         'Dnase I': [-0.30028264, 0.06592724, -1.6872665, -0.07824988, -0.36516234,
                           0.47394848, -0.16547703, 0.30670302, -0.32623452, -0.70181591,
                           -1.6872665, -0.75299878, -0.64558683, 0.03060385, -0.16547703,
                           1.10472336, 0.03925447, 0.10918037, -0.75299878, 0.53594464,
                           0.1928031, 2.27399978, 1.10472336, 2.14135683, -0.07392456,
                           -0.15250109, -0.07824988, 0.53594464, -1.1119998, -0.49059643,
                           0.30670302, 2.14135683, 0.55829209, -0.83445885, -0.32623452,
                           -0.07392456, -0.76164941, 0.16613034, -0.64558683, -1.1119998,
                           0.55829209, 0.91729311, -0.30028264, 0.03925447, -0.76164941,
                           0.06160193, -0.36516234, 0.1928031, 0.91729311, 4.52172102,
                           -0.70181591, -0.15250109, 0.06160193, 0.20577904, 0.03060385,
                           -0.49059643, -0.83445885, -1.10334918, 0.06592724, 0.10918037,
                           0.16613034, 0.20577904, 0.47394848, -2.61504624],
         'Dnase I-Rigid': [-0.08423873, -0.27526305, -1.02937647, 0.10745118, -0.42701756,
                                 0.1460554, -0.0463001, 0.72711543, 0.45488914, -0.76713403,
                                 -1.02937647, -0.31320167, -0.50289481, -0.3511403, -0.0463001,
                                 1.51583953, 0.29980668, -0.38907893, -0.31320167, 1.35676353,
                                 0.49349335, 2.11753285, 1.51583953, 2.15813384, 0.10745118,
                                 -0.12284295, 0.10745118, 1.35676353, -1.21507607, -1.58514408,
                                 0.72711543, 2.15813384, 0.96206868, -0.12284295, 0.45488914,
                                 0.10745118, -1.54853663, -0.27526305, -0.50289481, -1.21507607,
                                 0.96206868, 1.8759237, -0.08423873, 0.29980668, -1.54853663,
                                 -0.87961873, -0.42701756, 0.49349335, 1.8759237, -1.69563202,
                                 -0.76713403, -0.12284295, -0.87961873, -0.69192236, -0.3511403,
                                 -1.58514408, -0.12284295, -1.69563202, -0.27526305, -0.38907893,
                                 -0.27526305, -0.69192236, 0.1460554, 2.11753285],
         'MW-Daltons': [1., 1., 1., 1., -1., -1., -1., -1., 1., 1., 1., 1., -1.,
                       -1., -1., -1., 1., 1., 1., 1., -1., -1., -1., -1., 1., 1.,
                        1., 1., -1., -1., -1., -1., 1., 1., 1., 1., -1., -1., -1.,
                       -1., 1., 1., 1., 1., -1., -1., -1., -1., 1., 1., 1., 1.,
                       -1., -1., -1., -1., 1., 1., 1., 1., -1., -1., -1., -1.],
         'MW-kg':      [1., 1., 1., 1., -1., -1., -1., -1., 1., 1., 1., 1., -1.,
                       -1., -1., -1., 1., 1., 1., 1., -1., -1., -1., -1., 1., 1.,
                        1., 1., -1., -1., -1., -1., 1., 1., 1., 1., -1., -1., -1.,
                       -1., 1., 1., 1., 1., -1., -1., -1., -1., 1., 1., 1., 1.,
                       -1., -1., -1., -1., 1., 1., 1., 1., -1., -1., -1., -1.],
         'Nucleosome': [0.55530614, -0.50701865, 2.48680575, 0.2655812, 0.2655812,
                              -0.89331857, 0.2655812, 0.21729371, 1.27961849, 0.55530614,
                              2.48680575, 1.27961849, 0.79674359, -0.55530614, 0.2655812,
                              -0.55530614, 0.2655812, -0.7484561, 1.27961849, 0.45873116,
                              0.16900622, -2.34194328, -0.55530614, -2.00393085, 0.2655812,
                              0.16900622, 0.2655812, 0.45873116, 0.89331857, -0.98989355,
                              0.21729371, -2.00393085, -0.07243124, 1.66591842, 1.27961849,
                              0.2655812, -0.31386869, -0.7484561, 0.79674359, 0.89331857,
                              -0.07243124, 0.2655812, 0.55530614, 0.2655812, -0.31386869,
                              -1.27961849, 0.2655812, 0.16900622, 0.2655812, 0.2655812,
                              0.55530614, 0.16900622, -1.27961849, -1.37619348, -0.55530614,
                              -0.98989355, 1.66591842, 0.2655812, -0.50701865, -0.7484561,
                              -0.7484561, -1.37619348, -0.89331857, -2.34194328],
         'Nucleosome-Rigid': [-0.5618277, 0.49918947, -2.43267904, -0.27452139, -0.27452139,
                                    0.88982524, -0.27452139, -0.22663701, -1.27169272, -0.5618277,
                                    -2.43267904, -1.27169272, -0.79956948, 0.54791393, -0.27452139,
                                    0.54791393, -0.27452139, 0.74281178, -1.27169272, -0.46605893,
                                    -0.17875262, 2.38600227, 0.54791393, 2.03232988, -0.27452139,
                                    -0.17875262, -0.27452139, -0.46605893, -0.89449817, 0.98811424,
                                    -0.22663701, 2.03232988, 0.06234946, -1.64636703, -1.27169272,
                                    -0.27452139, 0.30429162, 0.74281178, -0.79956948, -0.89449817,
                                    0.06234946, -0.27452139, -0.5618277, -0.27452139, 0.30429162,
                                    1.2846614, -0.27452139, -0.17875262, -0.27452139, -0.27452139,
                                    -0.5618277, -0.17875262, 1.2846614, 1.38379048, 0.54791393,
                                    0.98811424, -1.64636703, -0.27452139, 0.49918947, 0.74281178,
                                    0.74281178, 1.38379048, 0.88982524, 2.38600227]
                     }
        self.properties = [
            "Bendability-DNAse",
            "Bendability-consensus",
            "Trinucleotide GC Content",
            "Nucleosome positioning",
            "Consensus_roll",
            "Consensus_Rigid",
            "Dnase I",
            "Dnase I-Rigid",
            "MW-Daltons",
            "MW-kg",
            "Nucleosome",
            "Nucleosome-Rigid",
        ]
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
            self.supInfo[props] = (originalValues - _mean) / _std

    def encode_into_kmers(self, seq):
        """
        Count kmer occurrences in a given seq.
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
        for l in range(len(seq) - 2):
            kmer = seq[l: l + 3].upper()
            try:
                # searches oligonucs for the specific kmer and return its index in that array
                index = self.oligonucs.index(kmer)
            except:
                print(
                    "The kmer {} is not in the possible oligonucs {}".format(
                        kmer, self.oligonucs
                    )
                )
            # append the kmer index number to create a list of kmers indexes (these indexes are in the oligonucs)
            seq_into_oligonucs_index.append(index)
            # add 1 to the count of what and how many kmers were found
            kmerCounts[index] += 1

        # find the normalized frequence for the kmers in the sequence
        normal_freq_of_kmers = kmerCounts / np.sum(kmerCounts)

        return seq_into_oligonucs_index, normal_freq_of_kmers

    def big_theta(self, index_kmer1, index_kmer2):
        ni = len(self.properties)  # ni = number of properties being extracted
        props_sum_for_kmers_at_i_and_ij = 0
        # sum of (p_v(kmer at pos i) - p_v(kmer at pos i+j))^2 for v in properties
        for v in self.properties:
            props_sum_for_kmers_at_i_and_ij += np.power(
                self.supInfo[v][index_kmer1] - self.supInfo[v][index_kmer2], 2
            )

        # divide the sum for ni = number of properties being extracted
        return props_sum_for_kmers_at_i_and_ij / ni

    def small_theta(self, seq_into_oligonucs_index):
        # [0, 1, 5, 22, 24, 35, 12, 51, 13]

        theta = np.zeros([self.lamb])

        for j in range(1, self.lamb + 1):
            # starts the usm with zero and adding the big_theta results
            big_theta_sum = 0
            for i in range(1, (self.L - j - 1)):
                # subtracting 1 to access the index of array seq_into_oligonucs_index
                big_theta_sum += self.big_theta(
                    seq_into_oligonucs_index[i - 1], seq_into_oligonucs_index[i + j - 1]
                )

            # j-1 is used to access since position 1 (equ. index=0) of the array
            theta[j - 1] = big_theta_sum / (self.L - j - 1)

        return theta

    def encode_into_pseknc(self, seq):
        seq_into_oligonucs_index, normal_freq_of_kmers = self.encode_into_kmers(seq)
        # w_theta is w * small_theta
        w_theta = self.w * self.small_theta(seq_into_oligonucs_index)
        feature_vector = np.concatenate((normal_freq_of_kmers, w_theta)) / (
                np.sum(normal_freq_of_kmers) + np.sum(w_theta)
        )
        return feature_vector

    def encode_fasta_into_pseknc(self, input_path, output_path, save=True):
        file_in = open(input_path, "r")
        if save:
            file_out = open(output_path, "w")
        else:
            encoded = []
        sigma = ""

        for i, line in enumerate(file_in):
            if line == "":
                continue
            elif line[0] == ">":
                sigma = line[1:-1]
            else:
                if line[-1:] == "\n":
                    line58 = line[:-1]
                else:
                    line58 = line

                if len(line58) != self.L:
                    print(
                        "Line {} has size {} instead of {}.".format(
                            i, len(line58), self.L
                        )
                    )
                    continue
                line_enconded = self.encode_into_pseknc(line58)
                if save:
                    file_out.write(str(list(line_enconded))[1:-1] + "," + self.synonyms(sigma) + "\n")
                else:
                    encoded.append(line_enconded)
        file_in.close()
        if save:
            file_out.close()
            return output_path
        else:
            return encoded
