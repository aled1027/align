#!/home/alex/anaconda3/bin/python
import unittest
import numpy as np
import align
from Bio import SeqIO

def _fasta_dict(filename):
    """Generates a dictionary mapping read name to sequence

    SeqIO.parse returns an iterator to records,
    where each record is a Bio.SeqIO.SeqRecord

    Args:
        filename (string): The path of fasta file.
    Returns
        Dictionary: Maps read names to their string
    """

    with open(filename, 'r') as fasta_file:
        ret_dict = {record.id: str(record.seq) \
                for record in SeqIO.parse(fasta_file, "fasta")}
    return ret_dict




class TestAlign(unittest.TestCase):
    def test_basic(self):
        """Verify that the aligner is putting out the correct score.
        Example on smith-waterman wikipedia page.
        """

        s0 = 'ACACACTA'
        s1 = 'AGCACACA'
        score, normalized_score, a0, a1 = align.align(s0, s1, local=True)

        self.assertEqual(score, 12)

    def test_align(self):
        d = _fasta_dict('example.fa')
        r0 = 'm150213_074729_42177R_c100777662550000001823160908051505_s1_p0/70715/9957_22166'
        r1 = 'm150126_093705_42156_c100779662550000001823165208251525_s1_p0/144605/28461_40297'
        size = 2000
        s0 = d[r0][:size]
        s1 = d[r1][:size]

        score, normalized_score, a0, a1 = align.align(s0, s1, local=False)
        print(normalized_score)

    def test_long_align(self):
        d = _fasta_dict('example.fa')
        r0 = 'm150213_074729_42177R_c100777662550000001823160908051505_s1_p0/70715/9957_22166'
        r1 = 'm150126_093705_42156_c100779662550000001823165208251525_s1_p0/144605/28461_40297'
        s0 = d[r0]
        s1 = d[r1]

        score, normalized_score, a0, a1 = align.align(s0, s1, local=False)
        print(score)
        print(normalized_score)


if __name__ == '__main__':
    unittest.main()




