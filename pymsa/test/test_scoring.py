import unittest

from pymsa.substitutionmatrix import PAM250, Blosum62
from pymsa.score import SumOfPairs, Star, Entropy, PercentageOfTotallyConservedColumns, PercentageOfNonGaps


class SumOfPairsTestCases(unittest.TestCase):
    def setUp(self):
        self.sumofpairs_PAM250 = SumOfPairs(PAM250())
        self.sumofpairs_Blosum62 = SumOfPairs(Blosum62())

    def tearDown(self):
        pass

    def test_basic_score_of_12_with_PAM250(self):
        # setup
        sequences = ['AA', 'AA', 'AA']

        # results
        result = self.sumofpairs_PAM250.compute(sequences)
        expected = 12

        # check
        self.assertEqual(expected, result)

    def test_basic_score_of_12_with_BLOSUM62(self):
        # setup
        sequences = ['AA', 'AA', 'AA']

        # results
        result = self.sumofpairs_Blosum62.compute(sequences)
        expected = 24

        # check
        self.assertEqual(expected, result)

    def test_basic_score_with_gaps_BLOSUM62(self):
        # setup
        sequences = ['FA', 'A-']

        # results
        result = self.sumofpairs_Blosum62.compute(sequences)
        expected = -10

        # check
        self.assertEqual(expected, result)

    def test_only_gaps_with_BLOSUM62(self):
        # setup
        sequences = ['---', '---']

        # results
        result = self.sumofpairs_Blosum62.compute(sequences)
        expected = 3

        # check
        self.assertEqual(expected, result)

    def test_get_score_of_column_with_only_gaps(self):
        # setup
        column = ['-', '-', '-']

        # results
        result = self.sumofpairs_Blosum62.get_score_of_k_column(column)
        expected = 3

        # check
        self.assertEqual(expected, result)

    def test_get_score_of_an_alignment(self):
        # setup
        sequences = \
            ['---GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYI------PPKGE----',
             '------MQDRVKRPMNAFIVWSRDQRRKMALENPRMR-NS-EISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRP---RRKAKMLPK',
             'MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMS-NL-DLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDH---PDLIQNAKK',
             '--------MHIKKPLNAFMLYMKEMRANVVAES-TLK-ESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK']
        # results
        result = self.sumofpairs_PAM250.compute(sequences)
        expected = 24

        # check
        self.assertEqual(expected, result)


class StarTestCases(unittest.TestCase):
    def setUp(self):
        self.star_PAM250 = Star(PAM250())
        self.star_Blosum62 = Star(Blosum62())

    def tearDown(self):
        pass

    def test_most_frequent_A_with_BLOSUM62(self):
        # setup
        sequences = ['AA', 'AC', 'AC']

        # results
        result = self.star_Blosum62.compute(sequences)
        expected = 30

        # check
        self.assertEqual(expected, result)

    def test_most_frequent_error_with_BLOSUM62(self):
        # setup
        sequences = ['AA', 'AC', 'AC']

        # results
        result = self.star_Blosum62.compute(sequences)
        expected = 22

        # check
        with self.assertRaises(Exception):
            self.assertEqual(expected, result)

    def test_most_frequent_with_PAM250(self):
        # setup
        sequences = ['AA', 'AC', 'AC']

        # results
        result = self.star_PAM250.compute(sequences)
        expected = 28

        # check
        self.assertEqual(expected, result)

    def test_most_frequent_error_with_PAM250(self):
        # setup
        sequences = ['AA', 'AC', 'AC']

        # results
        result = self.star_PAM250.compute(sequences)
        expected = 22

        # check
        with self.assertRaises(Exception):
            self.assertEqual(expected, result)

    def test_most_frequent_gaps_with_PAM250(self):
        # setup
        sequences = ['AA', 'A-', 'AC']

        # results
        result = self.star_PAM250.compute(sequences)
        expected = -2

        # check
        self.assertEqual(expected, result)

    def test_most_frequent_gaps_with_BLOSUM62(self):
        # setup
        sequences = ['AA', 'A-', 'AC']

        # results
        result = self.star_Blosum62.compute(sequences)
        expected = 8

        # check
        self.assertEqual(expected, result)


class EntropyTestCases(unittest.TestCase):
    def setUp(self):
        self.ent = Entropy()

    def tearDown(self):
        pass

    def test_get_entropy_of_a_column_with_gaps(self):
        # setup
        d = {"-": 0.8, "A": 0.2}

        # results
        result = round(self.ent.get_column_entropy(d), 1)
        expected = -0.5

        # check
        self.assertEqual(expected, result)

    def test_get_entropy_of_a_gapped_column(self):
        # setup
        d = {"-": 1}

        # results
        expected = 0
        result = self.ent.get_column_entropy(d)

        # check
        self.assertEqual(expected, result)

    def test_get_dictionary_of_a_five_letter_column(self):
        # setup
        column = ["-", "A", "C", "G", "T"]
        tot_seqs = len(column)

        # results
        expected = {"-": 1 / tot_seqs, "A": 1 / tot_seqs, "C": 1 / tot_seqs, "G": 1 / tot_seqs, "T": 1 / tot_seqs}
        result = self.ent.get_words_frequencies(column)

        # check
        self.assertEqual(expected, result)

    def test_get_dictionary_of_a_gapped_column(self):
        # setup
        column = ["-", "-", "-", "-", "-"]

        # results
        expected = {"-": 1}
        result = self.ent.get_words_frequencies(column)

        # check
        self.assertEqual(expected, result)

    def test_compute_of_four_seqs_with_no_gaps(self):
        # setup
        sequences = ["ACGT", "ACGT", "TGCA", "TGCA"]

        # results
        expected = -2.77
        result = round(self.ent.compute(sequences), 2)

        # check
        self.assertEqual(expected, result)

    def test_compute_of_three_seqs_with_gaps(self):
        # setup
        sequences = ["A-TGCAAT-G", "-CT-CCAT-A", "-TTAT-CTG-"]

        # results
        expected = -6.94
        result = round(self.ent.compute(sequences), 2)

        # check
        self.assertEqual(expected, result)

    def test_compute_of_two_gapped_seqs(self):
        # setup
        sequences = ["-----", "-----"]

        # results
        expected = 0
        result = self.ent.compute(sequences)

        # check
        self.assertEqual(expected, result)


class PercentageOfTotallyConservedColumnsTestCases(unittest.TestCase):
    def setUp(self):
        self.per = PercentageOfTotallyConservedColumns()

    def tearDown(self):
        pass

    def test_percentage_of_totally_conserved_columns_100(self):
        # setup
        sequences = ["AdddAAA", "AdddAAA"]

        # results
        result = self.per.compute(sequences)
        expected = 100.0

        # check
        self.assertEqual(result, expected)

    def test_percentage_of_totally_conserved_columns_50(self):
        # setup
        sequences = ["AB", "AC", "AC"]

        # results
        result = self.per.compute(sequences)
        expected = 50.0

        # check
        self.assertEqual(result, expected)

    def test_percentage_of_totally_conserved_columns_0(self):
        # setup
        sequences = ["ABCD", "DCBA"]

        # results
        result = self.per.compute(sequences)
        expected = 0.0

        # check
        self.assertEqual(result, expected)


class PercentageOfNonGapsTestCases(unittest.TestCase):
    def setUp(self):
        self.per = PercentageOfNonGaps()

    def tearDown(self):
        pass

    def test_percentage_of_non_gaps_100(self):
        # setup
        sequences = ["AB", "AC", "AC"]

        # results
        result = self.per.compute(sequences)
        expected = 100.0

        # check
        self.assertEqual(result, expected)

    def test_percentage_of_non_gaps_50(self):
        # setup
        sequences = ["A-", "A-"]

        # results
        result = self.per.compute(sequences)
        expected = 50.0

        # check
        self.assertEqual(result, expected)

    def test_percentage_of_non_gaps_0(self):
        # setup
        sequences = ["----", "----"]

        # results
        result = self.per.compute(sequences)
        expected = 0.0

        # check
        self.assertEqual(result, expected)

if __name__ == "__main__":
    unittest.main()