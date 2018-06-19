import unittest
import os

from pymsa.core.substitution_matrix import SubstitutionMatrix, FileMatrix, PAM250, Blosum62


class SubstitutionMatrixTestCases(unittest.TestCase):

    def test_should_default_gap_penalty_be_minus_eight(self):
        matrix = SubstitutionMatrix(-8, '-')

        self.assertEqual(-8, matrix.gap_penalty)

    def test_should_constructor__modify_the_gap_penalty(self):
        matrix = SubstitutionMatrix(-10, '-')

        self.assertEqual(-10, matrix.gap_penalty)

    def test_should_get_distance_return_the_gap_penalty_if_a_char_is_a_gap(self):
        matrix = SubstitutionMatrix(-8, '-')

        self.assertEqual(matrix.gap_penalty, matrix.get_distance('A', '-'))
        self.assertEqual(matrix.gap_penalty, matrix.get_distance('-', 'B'))

    def test_should_get_distance_return_one_if_the_two_chars_are_gaps(self):
        matrix = SubstitutionMatrix(-8, '-')

        self.assertEqual(1, matrix.get_distance('-', '-'))


class FileMatrixTestCases(unittest.TestCase):

    def test_should_default_gap_penalty_be_minus_eight(self):
        matrix = FileMatrix(path_to_file=os.path.dirname(__file__)+'/test_matrix.txt')

        self.assertEqual(matrix.gap_penalty, matrix.get_distance('A', '-'))
        self.assertEqual(matrix.gap_penalty, matrix.get_distance('-', 'B'))

    def test_should_get_distance_return_the_correct_values_if_there_are_no_gaps(self):
        matrix = FileMatrix(path_to_file=os.path.dirname(__file__)+'/test_matrix.txt')

        self.assertEqual(-1, matrix.get_distance('A', 'R'))
        self.assertEqual(-1, matrix.get_distance('R', 'A'))
        self.assertEqual(-3, matrix.get_distance('X', 'C'))
        self.assertEqual(+4, matrix.get_distance('I', 'I'))
        self.assertEqual(+4, matrix.get_distance('V', 'V'))


class PAM250TestCases(unittest.TestCase):

    def test_should_default_gap_penalty_be_minus_eight(self):
        matrix = PAM250()

        self.assertEqual(-8, matrix.gap_penalty)

    def test_should_constructor__modify_the_gap_penalty(self):
        matrix = PAM250(-10)

        self.assertEqual(-10, matrix.gap_penalty)

    def test_should_get_distance_return_the_gap_penalty_if_a_char_is_a_gap(self):
        matrix = PAM250()

        self.assertEqual(matrix.gap_penalty, matrix.get_distance('A', '-'))
        self.assertEqual(matrix.gap_penalty, matrix.get_distance('-', 'B'))

    def test_should_get_distance_return_one_if_the_two_chars_are_gaps(self):
        matrix = PAM250()

        self.assertEqual(1, matrix.get_distance('-', '-'))

    def test_should_get_distance_return_the_correct_values_if_there_are_no_gaps(self):
        matrix = PAM250()

        self.assertEqual(-2, matrix.get_distance('A', 'R'))
        self.assertEqual(-3, matrix.get_distance('N', 'F'))
        self.assertEqual(-3, matrix.get_distance('X', 'C'))
        self.assertEqual(+5, matrix.get_distance('I', 'I'))
        self.assertEqual(+4, matrix.get_distance('V', 'V'))

    def test_should_get_distance_throw_an_exception_if_a_char_is_invalid(self):
        matrix = PAM250()

        with self.assertRaises(Exception):
            matrix.get_distance('J', 'A')


class Blosum62TestCases(unittest.TestCase):

    def test_should_default_gap_penalty_be_minus_eight(self):
        matrix = Blosum62()

        self.assertEqual(-8, matrix.gap_penalty)

    def test_should_constructor__modify_the_gap_penalty(self):
        matrix = Blosum62(-10)

        self.assertEqual(-10, matrix.gap_penalty)

    def test_should_get_distance_return_the_gap_penalty_if_a_char_is_a_gap(self):
        matrix = Blosum62()

        self.assertEqual(matrix.gap_penalty, matrix.get_distance('A', '-'))
        self.assertEqual(matrix.gap_penalty, matrix.get_distance('-', 'B'))

    def test_should_get_distance_return_one_if_the_two_chars_are_gaps(self):
        matrix = Blosum62()

        self.assertEqual(1, matrix.get_distance('-', '-'))

    def test_should_get_distance_return_the_correct_values_if_there_are_no_gaps(self):
        matrix = Blosum62()

        self.assertEqual(-1, matrix.get_distance('A', 'R'))
        self.assertEqual(-3, matrix.get_distance('N', 'F'))
        self.assertEqual(-2, matrix.get_distance('X', 'C'))
        self.assertEqual(+4, matrix.get_distance('I', 'I'))
        self.assertEqual(+4, matrix.get_distance('V', 'V'))

    def test_should_get_distance_throw_an_exception_if_a_char_is_invalid(self):
        matrix = Blosum62()

        with self.assertRaises(Exception):
            matrix.get_distance('J', 'A')


if __name__ == '__main__':
    unittest.main()
