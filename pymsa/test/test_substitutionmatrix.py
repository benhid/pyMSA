import unittest

from pymsa.substitutionmatrix import SubstitutionMatrix, PAM250, Blosum62


class SubstitutionMatrixTestCases(unittest.TestCase):
    def setUp(self):
        print("setUp: RUNNING SubstutionMatrixTestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    def test_should_default_gap_penalty_be_minus_eight(self):
        matrix = SubstitutionMatrix()

        self.assertEqual(-8, matrix.get_gap_penalty())

    def test_should_constructor__modify_the_gap_penalty(self):
        matrix = SubstitutionMatrix(-10)

        self.assertEqual(-10, matrix.get_gap_penalty())

    def test_should_get_distance_return_the_gap_penalty_if_a_char_is_a_gap(self):
        matrix = SubstitutionMatrix()

        self.assertEqual(matrix.get_gap_penalty(), matrix.get_distance('A', '-'))
        self.assertEqual(matrix.get_gap_penalty(), matrix.get_distance('-', 'B'))

    def test_should_get_distance_return_one_if_the_two_chars_are_gaps(self):
        matrix = SubstitutionMatrix()

        self.assertEqual(1, matrix.get_distance('-', '-'))


class PAM250TestCases(unittest.TestCase):
    def setUp(self):
        print("setUp: RUNNING PAM250TestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    def test_should_default_gap_penalty_be_minus_eight(self):
        matrix = PAM250()

        self.assertEqual(-8, matrix.get_gap_penalty())

    def test_should_constructor__modify_the_gap_penalty(self):
        matrix = PAM250(-10)

        self.assertEqual(-10, matrix.get_gap_penalty())

    def test_should_get_distance_return_the_gap_penalty_if_a_char_is_a_gap(self):
        matrix = PAM250()

        self.assertEqual(matrix.get_gap_penalty(), matrix.get_distance('A', '-'))
        self.assertEqual(matrix.get_gap_penalty(), matrix.get_distance('-', 'B'))

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
    def setUp(self):
        print("setUp: RUNNING Blosum62TestCases")

    def tearDown(self):
        print("tearDown: TEST ENDED")

    def test_should_default_gap_penalty_be_minus_eight(self):
        matrix = Blosum62()

        self.assertEqual(-8, matrix.get_gap_penalty())

    def test_should_constructor__modify_the_gap_penalty(self):
        matrix = Blosum62(-10)

        self.assertEqual(-10, matrix.get_gap_penalty())

    def test_should_get_distance_return_the_gap_penalty_if_a_char_is_a_gap(self):
        matrix = Blosum62()

        self.assertEqual(matrix.get_gap_penalty(), matrix.get_distance('A', '-'))
        self.assertEqual(matrix.get_gap_penalty(), matrix.get_distance('-', 'B'))

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
