import unittest

from pymsa.core.msa import MSA
from pymsa.core.substitution_matrix import PAM250, Blosum62
from pymsa.core.score import Score, SumOfPairs, Star, Entropy, PercentageOfTotallyConservedColumns, PercentageOfNonGaps


class ScoreTestCases(unittest.TestCase):

    def test_should_raise_exception_if_sequences_have_different_lengths(self):
        # setup
        sequences = MSA(['AA', 'A', 'AA'])

        # check
        with self.assertRaises(Exception):
            Score(sequences)

    def test_should_raise_exception_if_msa_contains_only_one_sequence(self):
        # setup
        sequence = MSA(['AA'])

        # check
        with self.assertRaises(Exception):
            Score(sequence)


class SumOfPairsTestCases(unittest.TestCase):

    def test_basic_score_of_12_with_PAM250(self):
        # setup
        sequences = MSA(['AA', 'AA', 'AA'])

        # results
        result = SumOfPairs(sequences, PAM250()).compute()
        expected = 12

        # check
        self.assertEqual(expected, result)

    def test_basic_score_of_12_with_BLOSUM62(self):
        # setup
        sequences = MSA(['AA', 'AA', 'AA'])

        # results
        result = SumOfPairs(sequences, Blosum62()).compute()
        expected = 24

        # check
        self.assertEqual(expected, result)

    def test_basic_score_with_gaps_BLOSUM62(self):
        # setup
        sequences = MSA(['FA', 'A-'])

        # results
        result = SumOfPairs(sequences, Blosum62()).compute()
        expected = -10

        # check
        self.assertEqual(expected, result)

    def test_only_gaps_with_BLOSUM62(self):
        # setup
        sequences = MSA(['---', '---'])

        # results
        result = SumOfPairs(sequences, Blosum62()).compute()
        expected = 3

        # check
        self.assertEqual(expected, result)

    def test_get_score_of_column_with_only_gaps(self):
        # setup
        column = MSA(['-', '-', '-'])

        # results
        result = SumOfPairs(column).get_column_score(0)
        expected = 3

        # check
        self.assertEqual(expected, result)

    def test_get_score_of_an_alignment(self):
        # setup
        sequences = \
            MSA(['---GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYI------PPKGE----',
             '------MQDRVKRPMNAFIVWSRDQRRKMALENPRMR-NS-EISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRP---RRKAKMLPK',
             'MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMS-NL-DLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDH---PDLIQNAKK',
             '--------MHIKKPLNAFMLYMKEMRANVVAES-TLK-ESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK'])

        # results
        result = SumOfPairs(sequences, PAM250()).compute()
        expected = 24

        # check
        self.assertEqual(expected, result)


class StarTestCases(unittest.TestCase):

    def test_most_frequent_A_with_BLOSUM62(self):
        # setup
        sequences = MSA(['AA', 'AC', 'AC'])

        # results
        result = Star(sequences, Blosum62()).compute()
        expected = 30

        # check
        self.assertEqual(expected, result)

    def test_most_frequent_with_PAM250(self):
        # setup
        sequences = MSA(['AA', 'AC', 'AC'])

        # results
        result = Star(sequences, PAM250()).compute()
        expected = 28

        # check
        self.assertEqual(expected, result)

    def test_most_frequent_gaps_with_PAM250(self):
        # setup
        sequences = MSA(['AA', 'A-', 'AC'])

        # results
        result = Star(sequences, PAM250()).compute()
        expected = -2

        # check
        self.assertEqual(expected, result)

    def test_most_frequent_gaps_with_BLOSUM62(self):
        # setup
        sequences = MSA(['AA', 'A-', 'AC'])

        # results
        result = Star(sequences, Blosum62()).compute()
        expected = 8

        # check
        self.assertEqual(expected, result)


class EntropyTestCases(unittest.TestCase):

    def test_get_entropy_of_a_column_with_gaps(self):
        # setup
        d = {"-": 0.8, "A": 0.2}

        # results
        result = round(Entropy.get_entropy(d), 1)
        expected = -0.5

        # check
        self.assertEqual(expected, result)

    def test_get_entropy_of_a_gapped_column(self):
        # setup
        d = {"-": 1}

        # results
        expected = 0
        result = Entropy.get_entropy(d)

        # check
        self.assertEqual(expected, result)

    def test_get_dictionary_of_a_five_letter_column(self):
        # setup
        column = ["-", "A", "C", "G", "T"]
        tot_seqs = len(column)

        # results
        expected = {"-": 1 / tot_seqs, "A": 1 / tot_seqs, "C": 1 / tot_seqs, "G": 1 / tot_seqs, "T": 1 / tot_seqs}
        result = Entropy.get_words_frequencies(column)

        # check
        self.assertEqual(expected, result)

    def test_get_dictionary_of_a_gapped_column(self):
        # setup
        column = ["-", "-", "-", "-", "-"]

        # results
        expected = {"-": 1}
        result = Entropy.get_words_frequencies(column)

        # check
        self.assertEqual(expected, result)

    def test_compute_of_four_seqs_with_no_gaps(self):
        # setup
        sequences = MSA(["ACGT", "ACGT", "TGCA", "TGCA"])

        # results
        expected = -2.77
        result = round(Entropy(sequences).compute(), 2)

        # check
        self.assertEqual(expected, result)

    def test_compute_of_three_seqs_with_gaps(self):
        # setup
        sequences = MSA(["A-TGCAAT-G", "-CT-CCAT-A", "-TTAT-CTG-"])

        # results
        expected = -6.94
        result = round(Entropy(sequences).compute(), 2)

        # check
        self.assertEqual(expected, result)

    def test_compute_of_two_gapped_seqs(self):
        # setup
        sequences = MSA(["-----", "-----"])

        # results
        expected = 0
        result = Entropy(sequences).compute()

        # check
        self.assertEqual(expected, result)


class PercentageOfTotallyConservedColumnsTestCases(unittest.TestCase):

    def test_percentage_of_totally_conserved_columns_100(self):
        # setup
        sequences = MSA(["AdddAAA", "AdddAAA"])

        # results
        result = PercentageOfTotallyConservedColumns(sequences).compute()
        expected = 100.0

        # check
        self.assertEqual(result, expected)

    def test_percentage_of_totally_conserved_columns_50(self):
        # setup
        sequences = MSA(["AB", "AC", "AC"])

        # results
        result = PercentageOfTotallyConservedColumns(sequences).compute()
        expected = 50.0

        # check
        self.assertEqual(result, expected)

    def test_percentage_of_totally_conserved_columns_0(self):
        # setup
        sequences = MSA(["ABCD", "DCBA"])

        # results
        result = PercentageOfTotallyConservedColumns(sequences).compute()
        expected = 0.0

        # check
        self.assertEqual(result, expected)


class PercentageOfNonGapsTestCases(unittest.TestCase):

    def test_percentage_of_non_gaps_100(self):
        # setup
        sequences = MSA(["AB", "AC", "AC"])

        # results
        result = PercentageOfNonGaps(sequences).compute()
        expected = 100.0

        # check
        self.assertEqual(result, expected)

    def test_percentage_of_non_gaps_50(self):
        # setup
        sequences = MSA(["A-", "A-"])

        # results
        result = PercentageOfNonGaps(sequences).compute()
        expected = 50.0

        # check
        self.assertEqual(result, expected)

    def test_percentage_of_non_gaps_0(self):
        # setup
        sequences = MSA(["----", "----"])

        # results
        result = PercentageOfNonGaps(sequences).compute()
        expected = 0.0

        # check
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
