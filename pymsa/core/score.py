import itertools
import subprocess
import logging
import math
import os
import urllib.request

from urllib.error import HTTPError
from collections import Counter
from pathlib import Path
import urllib.request

from pymsa.core.substitutionmatrix import SubstitutionMatrix, PAM250

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Score:

    def __init__(self):
        pass

    def compute(self, align_sequences: list) -> float:
        """ Compute the core.
        :param align_sequences: List of sequences (as strings).
        :return: Value of the core.
        """
        if not all(len(sequence) == len(align_sequences[0]) for sequence in align_sequences):
            raise Exception("All the sequences in the FASTA file must be aligned!")

        return self._compute(align_sequences)

    def _compute(self, sequences: list) -> float:
        pass

    def _get_score_of_two_chars(self, substitution_matrix: SubstitutionMatrix, char_a: str, char_b: str) -> int:
        """ Return the core of two chars using the substituion matrix.
        :param substitution_matrix: Matrix of scores such as PAM250, Blosum62, etc.
        :param char_a: First char.
        :param char_b: Second char.
        :return: Value of the core.
        """
        return int(substitution_matrix.get_distance(char_a, char_b))

    def is_minimization(self) -> bool:
        pass

    def get_name(self) -> str:
        return type(self).__name__


class Entropy(Score):

    def _compute(self, align_sequences: list) -> float:
        length_of_sequence = len(align_sequences[0])
        column = []
        final_score = 0

        for k in range(length_of_sequence):
            for sequence in align_sequences:
                column.append(sequence[k])
            column_chars_and_frequencies = self.get_words_frequencies(column)
            final_score += self.get_column_entropy(column_chars_and_frequencies)
            column.clear()

        return final_score

    def get_words_frequencies(self, words: list) -> dict:
        """ Get dictionary of words frequencies of a list of words. """
        word_freq = [words.count(w) / len(words) for w in words]
        return dict(zip(words, word_freq))

    def get_column_entropy(self, column: dict) -> float:
        """ Calculates the Minimum Entropy for the current column. """
        current_entropy = 0

        for key, value in column.items():
            current_entropy += value * math.log(value)
            logger.debug('Character {0}, frecuency: {1}, entropy: {2})'.format(key, value, value * math.log(value)))

        return current_entropy

    def is_minimization(self) -> bool:
        return False


class Star(Score):

    def __init__(self, substitution_matrix: SubstitutionMatrix = PAM250()):
        super().__init__()
        self.substitution_matrix = substitution_matrix

    def _compute(self, align_sequences: list) -> int:
        length_of_sequence = len(align_sequences[0])
        column = []
        final_score = 0

        for k in range(length_of_sequence):
            for sequence in align_sequences:
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
            score_of_column += self._get_score_of_two_chars(self.substitution_matrix, most_frequent_char, char)

        logger.debug('Score of column: {0}'.format(score_of_column))
        return score_of_column

    def is_minimization(self) -> bool:
        return False


class SumOfPairs(Score):

    def __init__(self, substitution_matrix: SubstitutionMatrix = PAM250()):
        super().__init__()
        self.substitution_matrix = substitution_matrix

    def _compute(self, align_sequences: list) -> int:
        length_of_sequence = len(align_sequences[0])
        column = []
        final_score = 0

        for k in range(length_of_sequence):
            for sequence in align_sequences:
                column.append(sequence[k])
            final_score += self.get_score_of_k_column(column)
            column.clear()

        return final_score

    def get_score_of_k_column(self, column: list) -> int:
        """ Compare the most frequent element of the list 'column' with the others only one time. """
        score_of_column = 0

        for char_a, char_b in self.possible_combinations(column):
            score_of_column += self._get_score_of_two_chars(self.substitution_matrix, char_a, char_b)

        logger.debug('Score of column: {0}'.format(score_of_column))
        return score_of_column

    def possible_combinations(self, column) -> itertools.combinations:
        return itertools.combinations(column, 2)

    def is_minimization(self) -> bool:
        return False


class PercentageOfNonGaps(Score):

    def _compute(self, align_sequences: list) -> float:
        length_of_sequence = len(align_sequences[0])
        no_of_gaps = 0

        for sequence in align_sequences:
            no_of_gaps += sequence.count('-')

        logger.debug('Total number of gaps: {0}'.format(no_of_gaps))
        logger.debug('Total number of non-gaps: {0}'.format(length_of_sequence * 2 - no_of_gaps))

        return 100 - (no_of_gaps / (length_of_sequence * len(align_sequences)) * 100)

    def is_minimization(self) -> bool:
        return True


class PercentageOfTotallyConservedColumns(Score):

    def _compute(self, align_sequences: list) -> float:
        length_sequence = len(align_sequences[0])
        no_of_conserved_columns = 0
        column = []

        for k in range(length_sequence):
            for sequence in align_sequences:
                column.append(sequence[k])
            if len(set(column)) <= 1:
                no_of_conserved_columns += 1
            column.clear()

        logger.debug(
            'Total number of conserved columns: {0} out of {1}'.format(no_of_conserved_columns, length_sequence))

        return no_of_conserved_columns / length_sequence * 100

    def is_minimization(self) -> bool:
        return False


class Strike:

    def __init__(self):
        self.out_connection_path = 'strike/in.con'
        self.out_alignment_path = 'strike/aln.fa'

    def compute(self, align_sequences: list, sequences_id: list, chains: list):
        return self._compute(align_sequences, sequences_id, chains)

    def _compute(self, align_sequences: list, sequences_id: list, chains: list) -> float:
        # Check if directory exists (otherwise, create it)
        os.makedirs('strike/', exist_ok=True)

        if not Path(self.out_connection_path).is_file() and not Path(self.out_alignment_path).is_file():
            with open(self.out_connection_path, "w+") as c_file, open(self.out_alignment_path, "w+") as a_file:
                for i in range(len(align_sequences)):
                    self.get_pdb(sequences_id[i])
                    c_file.writelines(sequences_id[i] + ' ./strike/' + sequences_id[i] + '.pdb ' + chains[i] + '\n')
                    a_file.writelines('>' + sequences_id[i] + '\n' + align_sequences[i] + '\n')

        bytes = subprocess.check_output('strike -c ' + self.out_connection_path + ' -a ' + self.out_alignment_path,
                                        shell=True, env=os.environ.copy())
        score = "".join(map(chr, bytes)).split('\n')[-2]

        return float(score)

    def get_pdb(self, pdb_id: str) -> str:
        base_url = 'https://files.rcsb.org/download/'
        base_url += pdb_id + '.pdb'

        out_pdb_path = 'strike/' + pdb_id + '.pdb'

        if not Path(out_pdb_path).is_file():
            logger.debug('Downloading .pdb file...')

            try:
                req = urllib.request.Request(base_url)
                f = urllib.request.urlopen(req)
            except HTTPError:
                raise Exception("PDB not found in rcsb")

            with open(out_pdb_path, 'w+') as pdb_file:
                pdb_file.write(f.read().decode('unicode_escape'))

        return out_pdb_path

    def is_minimization(self) -> bool:
        return False
