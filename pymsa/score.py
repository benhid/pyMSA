import itertools
import logging
import math
from collections import Counter

from pymsa.substitutionmatrix import SubstitutionMatrix, PAM250

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Score:
    """ Class representing MSA (Multiple Sequence Alignment) scores
    
    Requirements:
    - All the sequences in an msa must be aligned
    - The gap character is '-'
    """

    def __init__(self):
        pass

    def compute(self, sequences: list) -> float:
        """ Compute the score 
        
        :param sequences: list of sequences (as strings)
        :return: the value of the score
        """
        pass

    def get_score_of_two_chars(self, substitution_matrix: SubstitutionMatrix, char_a: str, char_b: str) -> int:
        """ Return the score of two chars using the substituion matrix.
        :param substitution_matrix: Matrix of scores such as PAM250, Blosum62, etc.
        :param charA: First char.
        :param charB: Second char.
        :return: Value of the score.
        """
        return int(substitution_matrix.get_distance(char_a, char_b))

    def get_name(self) -> str:
        return type(self).__name__

    def _raiser(self, e): raise Exception(e)


class Entropy(Score):

    def compute(self, sequences: list) -> float:
        length_of_sequence = len(sequences[0])
        column = []
        final_score = 0

        for k in range(length_of_sequence):
            for sequence in sequences:
                column.append(sequence[k])
            column_chars_and_frequencies = self.get_words_frequencies(column)
            final_score += self.get_column_entropy(column_chars_and_frequencies)
            column.clear()

        return final_score

    def get_words_frequencies(self, words: list) -> dict:
        """ Get dictionary of words frequencies of a list of words. """
        word_freq = [words.count(w)/len(words) for w in words]
        return dict(zip(words, word_freq))

    def get_column_entropy(self, column: dict) -> float:
        """Calculates the Minimum Entropy for the current column. """
        current_entropy = 0

        for key, value in column.items():
            current_entropy += value * math.log(value)
            logger.debug('Character {0}, frecuency: {1}, entropy: {2})'.format(key, value, value * math.log(value)))

        return current_entropy


class Star(Score):

    def __init__(self, substitution_matrix: SubstitutionMatrix = PAM250()):
        super().__init__()
        self.substitution_matrix = substitution_matrix

    def compute(self, sequences: list) -> int:
        length_of_sequence = len(sequences[0])
        column = []
        final_score = 0

        for k in range(length_of_sequence):
            for sequence in sequences:
                column.append(sequence[k])
            final_score += self.get_score_of_k_column(column)
            column.clear()

        return final_score

    def get_score_of_k_column(self, column: list) -> int:
        """ Compare the most frequent element of the list 'column' with the others only one time.

        :param column: List of chars.
        :return: Score of two chars.
        """

        score_of_column = 0
        most_frequent_char = Counter(column).most_common(1)[0][0]

        for char in column:
            score_of_column += self.get_score_of_two_chars(self.substitution_matrix, most_frequent_char, char)

        logger.debug('Score of column: {0}'.format(score_of_column))
        return score_of_column


class SumOfPairs(Score):

    def __init__(self, substitution_matrix: SubstitutionMatrix = PAM250()):
        super().__init__()
        self.substitution_matrix = substitution_matrix

    def compute(self, sequences: list) -> int:
        length_of_sequence = len(sequences[0])
        column = []
        final_score = 0

        for k in range(length_of_sequence):
            for sequence in sequences:
                column.append(sequence[k])
            final_score += self.get_score_of_k_column(column)
            column.clear()

        return final_score

    def get_score_of_k_column(self, column: list) -> int:
        """ Compare the most frequent element of the list 'column' with the others only one time. """
        score_of_column = 0

        for char_a, char_b in self.possible_combinations(column):
            score_of_column += self.get_score_of_two_chars(self.substitution_matrix, char_a, char_b)

        logger.debug('Score of column: {0}'.format(score_of_column))
        return score_of_column

    def possible_combinations(self, column):
        return itertools.combinations(column, 2)


class PercentageOfNonGaps(Score):

    def compute(self, sequences: list) -> float:
        length_of_sequence = len(sequences[0])
        no_of_gaps = 0

        for sequence in sequences:
            no_of_gaps += sequence.count('-')

        logger.debug('Total number of gaps: {0}'.format(no_of_gaps))
        logger.debug('Total number of non-gaps: {0}'.format(length_of_sequence*2 - no_of_gaps))

        return 100 - (no_of_gaps / (length_of_sequence * len(sequences)) * 100)


class PercentageOfTotallyConservedColumns(Score):

    def compute(self, sequences: list) -> float:
        length_sequence = len(sequences[0])
        no_of_conserved_columns = 0
        column = []

        for k in range(length_sequence):
            for sequence in sequences:
                column.append(sequence[k])
            if len(set(column)) <= 1:
                no_of_conserved_columns += 1
            column.clear()

        logger.debug('Total number of conserved columns: {0} out of {1}'.format(no_of_conserved_columns, length_sequence))

        return no_of_conserved_columns / length_sequence * 100
