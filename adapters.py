"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This module contains the class and sequences for known adapters used in Oxford Nanopore library
preparation kits.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""


class Adapter(object):

    def __init__(self, name, start_sequence=None, end_sequence=None, both_ends_sequence=None):
        self.name = name
        self.start_sequence = start_sequence if start_sequence else []
        self.end_sequence = end_sequence if end_sequence else []
        if both_ends_sequence:
            self.start_sequence = both_ends_sequence
            self.end_sequence = both_ends_sequence
        self.best_start_score, self.best_end_score = 0.0, 0.0

    def best_start_or_end_score(self):
        return max(self.best_start_score, self.best_end_score)

    def is_barcode(self):
        return self.name.startswith('Barcode ')

    def barcode_direction(self):
        if '_rev' in self.start_sequence[0]:
            return 'reverse'
        else:
            return 'forward'

    def get_barcode_name(self):
        """
        Gets the barcode name for the output files. We want a concise name, so it looks at all
        options and chooses the shortest.
        """
        possible_names = [self.name]
        if self.start_sequence:
            possible_names.append(self.start_sequence[0])
        if self.end_sequence:
            possible_names.append(self.end_sequence[0])
        barcode_name = sorted(possible_names, key=lambda x: len(x))[0]
        return barcode_name.replace(' ', '_')


# INSTRUCTIONS FOR ADDING CUSTOM ADAPTERS
# ---------------------------------------
# If you need Porechop to remove adapters that aren't included, you can add your own my modifying
# the ADAPTERS list below.
#
# Here is the format for a normal adapter:
#     Adapter('Adapter_set_name',
#             start_sequence=('Start_adapter_name', 'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'),
#             end_sequence=('End_adapter_name', 'AACCGGTTAACCGGTTAACCGGTTAACCGGTT'))
#
# You can exclude start_sequence and end_sequence as appropriate.
#
# If you have custom Barcodes, make sure that the adapter set name starts with 'Barcode '. Also,
# remove the existing barcode sequences from this file to avoid conflicts:
#     Adapter('Barcode 1',
#             start_sequence=('Barcode_1_start', 'AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT'),
#             end_sequence=('Barcode_1_end', 'AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT')),
#     Adapter('Barcode 2',
#             start_sequence=('Barcode_2_start', 'TTTTTTTTGGGGGGGGCCCCCCCCAAAAAAAA'),
#             end_sequence=('Barcode_2_end', 'TTTTTTTTGGGGGGGGCCCCCCCCAAAAAAAA'))


ADAPTERS = [Adapter('SQK-NSK007',
                    start_sequence=('SQK-NSK007_Y_Top', 'AATGTACTTCGTTCAGTTACGTATTGCT'),
                    end_sequence=('SQK-NSK007_Y_Bottom', 'GCAATACGTAACTGAACGAAGT')),


            Adapter('Rapid',
                    start_sequence=('Rapid_adapter',
                                    'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA')),

            Adapter('RBK004_upstream',
                    start_sequence=('RBK004_upstream', 'AATGTACTTCGTTCAGTTACGGCTTGGGTGTTTAACC')),


            Adapter('SQK-MAP006',
                    start_sequence=('SQK-MAP006_Y_Top_SK63',    'GGTTGTTTCTGTTGGTGCTGATATTGCT'),
                    end_sequence=  ('SQK-MAP006_Y_Bottom_SK64', 'GCAATATCAGCACCAACAGAAA')),

            Adapter('SQK-MAP006 short',
                    start_sequence=('SQK-MAP006_Short_Y_Top_LI32',    'CGGCGTCTGCTTGGGTGTTTAACCT'),
                    end_sequence=  ('SQK-MAP006_Short_Y_Bottom_LI33', 'GGTTAAACACCCAAGCAGACGCCG')),


            # The PCR adapters are used both in PCR DNA kits and some cDNA kits.
            Adapter('PCR adapters 1',
                    start_sequence=('PCR_1_start', 'ACTTGCCTGTCGCTCTATCTTC'),
                    end_sequence=  ('PCR_1_end',   'GAAGATAGAGCGACAGGCAAGT')),

            Adapter('PCR adapters 2',
                    start_sequence=('PCR_2_start', 'TTTCTGTTGGTGCTGATATTGC'),
                    end_sequence=  ('PCR_2_end',   'GCAATATCAGCACCAACAGAAA')),

            Adapter('PCR adapters 3',
                    start_sequence=('PCR_3_start', 'TACTTGCCTGTCGCTCTATCTTC'),
                    end_sequence=  ('PCR_3_end',   'GAAGATAGAGCGACAGGCAAGTA')),


            # 1D^2 kit adapters are interesting. ONT provided the following sequences on their site:
            #   start: GGCGTCTGCTTGGGTGTTTAACCTTTTTGTCAGAGAGGTTCCAAGTCAGAGAGGTTCCT
            #   end:   GGAACCTCTCTCTGACTTGGAACCTCTCTGACAAAAAGGTTAAACACCCAAGCAGACGCCAGCAAT
            # But when looking at actual reads, I found two parts. The first corresponds to one end
            # of the provided sequences (through slightly different):
            Adapter('1D^2 part 1',
                    start_sequence=('1D2_part_1_start', 'GAGAGGTTCCAAGTCAGAGAGGTTCCT'),
                    end_sequence=  ('1D2_part_1_end',   'AGGAACCTCTCTGACTTGGAACCTCTC')),
            # and the second part corresponds to the other end, combined with a bit of standard 1D
            # adapter:
            Adapter('1D^2 part 2',
                    start_sequence=('1D2_part_2_start', 'CTTCGTTCAGTTACGTATTGCTGGCGTCTGCTT'),
                    end_sequence=  ('1D2_part_2_end',   'CACCCAAGCAGACGCCAGCAATACGTAACT')),
            # The middle part of the provided sequences is less common, so I've left it out of the
            # adapter sequences here.


            Adapter('cDNA SSP',
                    start_sequence=('cDNA_SSP',     'TTTCTGTTGGTGCTGATATTGCTGCCATTACGGCCGGG'),
                    end_sequence=  ('cDNA_SSP_rev', 'CCCGGCCGTAATGGCAGCAATATCAGCACCAACAGAAA')),


            # Some barcoding kits (like the native barcodes) use the rev comp barcode at the start
            # of the read and the forward barcode at the end of the read.
           

            # Other barcoding kits (like the PCR and rapid barcodes) use the forward barcode at the
            # start of the read and the rev comp barcode at the end of the read.
            Adapter('Barcode 1 (forward)',
                    start_sequence=('F01', 'AAGAAAGTTGTCGGTGTCTTTGTG'),
                    end_sequence=('F01_rev', 'CACAAAGACACCGACAACTTTCTT')),
            Adapter('Barcode 2 (forward)',
                    start_sequence=('F02', 'TCGATTCCGTTTGTAGTCGTCTGT'),
                    end_sequence=('F02_rev', 'ACAGACGACTACAAACGGAATCGA')),
            Adapter('Barcode 4 (forward)',
                    start_sequence=('F04', 'TTCGGATTCTATCGTGTTTCCCTA'),
                    end_sequence=('F04_rev', 'TAGGGAAACACGATAGAATCCGAA')),
            Adapter('Barcode 5 (forward)',
                    start_sequence=('F05', 'CTTGTCCAGGGTTTGTGTAACCTT'),
                    end_sequence=('F05_rev', 'AAGGTTACACAAACCCTGGACAAG')),
            Adapter('Barcode 6 (forward)',
                    start_sequence=('F06', 'TTCTCGCAAAGGCAGAAAGTAGTC'),
                    end_sequence=('F06_rev', 'GACTACTTTCTGCCTTTGCGAGAA')),
            Adapter('Barcode 7 (forward)',
                    start_sequence=('F07', 'GTGTTACCGTGGGAATGAATCCTT'),
                    end_sequence=('F07_rev', 'AAGGATTCATTCCCACGGTAACAC')),
            Adapter('Barcode 8 (forward)',
                    start_sequence=('F08', 'TTCAGGGAACAAACCAAGTTACGT'),
                    end_sequence=('F08_rev', 'ACGTAACTTGGTTTGTTCCCTGAA')),
            Adapter('Barcode 9 (forward)',
                    start_sequence=('F09', 'AACTAGGCACAGCGAGTCTTGGTT'),
                    end_sequence=('F09_rev', 'AACCAAGACTCGCTGTGCCTAGTT'))]


def make_full_native_barcode_adapter(barcode_num):
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (reverse)'][0]
    start_barcode_seq = barcode.start_sequence[1]
    end_barcode_seq = barcode.end_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA' + start_barcode_seq + 'CAGCACCT'
    end_full_seq = 'AGGTGCTG' + end_barcode_seq + 'TTAACCTTAGCAATACGTAACTGAACGAAGT'

    return Adapter('Native barcoding ' + str(barcode_num) + ' (full sequence)',
                   start_sequence=('NB' + '%02d' % barcode_num + '_start', start_full_seq),
                   end_sequence=('NB' + '%02d' % barcode_num + '_end', end_full_seq))


def make_old_full_rapid_barcode_adapter(barcode_num):  # applies to SQK-RBK001
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (forward)'][0]
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACG' + 'TATTGCT' + start_barcode_seq + \
                     'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter('Rapid barcoding ' + str(barcode_num) + ' (full sequence, old)',
                   start_sequence=('RB' + '%02d' % barcode_num + '_full', start_full_seq))


def make_new_full_rapid_barcode_adapter(barcode_num):  # applies to SQK-RBK004
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (forward)'][0]
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACG' + 'GCTTGGGTGTTTAACC' + start_barcode_seq + \
                     'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter('Rapid barcoding ' + str(barcode_num) + ' (full sequence, new)',
                   start_sequence=('RB' + '%02d' % barcode_num + '_full', start_full_seq))
