import itertools
import logging
import math
import os
import urllib.request
from abc import abstractmethod, ABC
from collections import Counter
from pathlib import Path
from urllib.error import HTTPError

from pymsa.core.msa import MSA
from pymsa.core.substitution_matrix import SubstitutionMatrix, PAM250
from pymsa.util.tool import StrikeEx

LOGGER = logging.getLogger('pyMSA')


def get_score_of_two_chars(substitution_matrix: SubstitutionMatrix, char_a: str, char_b: str) -> int:
    """
    Return the core of two chars using the substitution matrix.

    :param substitution_matrix: Matrix of scores such as PAM250, Blosum62, etc.
    :param char_a: First char.
    :param char_b: Second char.
    :return: Value of the core.
    """
    return int(substitution_matrix.get_distance(char_a, char_b))


class Score(ABC):

    def __init__(self, msa: MSA):
        self.msa = msa
        assert self.msa.is_valid, 'MSA is not valid'

    def compute(self) -> float:
        final_score = 0

        for k in range(len(self.msa)):
            final_score += self.get_column_score(k)

        return final_score

    def get_column(self, k: int) -> list:
        return [seq[k] for seq in self.msa.sequences]

    @abstractmethod
    def get_column_score(self, k: int) -> float:
        pass

    @staticmethod
    @abstractmethod
    def is_minimization() -> bool:
        pass


class Entropy(Score):

    def get_column_score(self, k) -> float:
        column = self.get_column(k)
        column_chars_and_frequencies = self.get_words_frequencies(column)

        return self.get_entropy(column_chars_and_frequencies)

    @staticmethod
    def get_words_frequencies(words: list) -> dict:
        """
        Compute frequencies of words inside a list.

        :param words: List of words.
        :return: Dictionary with computed frequencies.
        """
        word_freq = [words.count(w) / len(words) for w in words]
        return dict(zip(words, word_freq))

    @staticmethod
    def get_entropy(column: dict) -> float:
        """
        Compute the Minimum Entropy for the current column.
        """
        current_entropy = 0

        for key, value in column.items():
            current_entropy += value * math.log(value)

        return current_entropy

    @staticmethod
    def is_minimization() -> bool:
        return False


class Star(Score):

    def __init__(self, msa: MSA, substitution_matrix: SubstitutionMatrix = PAM250()):
        super(Star, self).__init__(msa=msa)
        self.substitution_matrix = substitution_matrix

    def get_column_score(self, k: int) -> float:
        column = self.get_column(k)

        score_of_column = 0
        most_frequent_char = Counter(column).most_common(1)[0][0]

        for char in column:
            score_of_column += get_score_of_two_chars(self.substitution_matrix, most_frequent_char, char)

        return score_of_column

    @staticmethod
    def is_minimization() -> bool:
        return False


class SumOfPairs(Score):

    def __init__(self, msa: MSA, substitution_matrix: SubstitutionMatrix = PAM250()):
        super(SumOfPairs, self).__init__(msa=msa)
        self.substitution_matrix = substitution_matrix

    def get_column_score(self, k: int) -> float:
        column = self.get_column(k)

        score_of_column = 0
        for char_a, char_b in self.possible_combinations(column):
            score_of_column += get_score_of_two_chars(self.substitution_matrix, char_a, char_b)

        return score_of_column

    @staticmethod
    def possible_combinations(column) -> itertools.combinations:
        return itertools.combinations(column, 2)

    @staticmethod
    def is_minimization() -> bool:
        return False


class PercentageOfNonGaps(Score):

    def compute(self) -> float:
        final_score = 0

        for k in range(len(self.msa)):
            final_score += self.get_column_score(k)

        return 100 - (final_score / (len(self.msa) * self.msa.number_of_sequences) * 100)

    def get_column_score(self, k: int) -> float:
        column = self.get_column(k)
        return column.count('-')

    @staticmethod
    def is_minimization() -> bool:
        return True


class PercentageOfTotallyConservedColumns(Score):

    def compute(self) -> float:
        final_score = 0

        for k in range(len(self.msa)):
            final_score += self.get_column_score(k)

        return final_score / len(self.msa) * 100

    def get_column_score(self, k: int) -> float:
        column = self.get_column(k)
        conserved_column = 0

        if len(set(column)) <= 1:
            conserved_column = 1

        return conserved_column

    @staticmethod
    def is_minimization() -> bool:
        return False


class Strike:

    def __init__(self, aligned_sequences: list, exe_path: str = '/usr/local/bin/strike'):
        self.aligned_sequences = aligned_sequences
        self.no_sequences = len(self.aligned_sequences)
        self.length = len(self.aligned_sequences[0])

        self.out_connection_path = os.path.abspath('strike/in.con')
        self.out_alignment_path = os.path.abspath('strike/aln.fa')
        self.exe_path = os.path.abspath(exe_path)

    def compute(self, sequences_id: list, chains: list) -> float:
        return self.evaluate(sequences_id, chains)

    def evaluate(self, sequences_id: list, chains: list) -> float:
        os.makedirs(os.path.abspath('strike'), exist_ok=True)

        if not Path(self.out_connection_path).is_file() and not Path(self.out_alignment_path).is_file():
            with open(self.out_connection_path, 'w+') as c_file, open(self.out_alignment_path, 'w+') as a_file:
                for i in range(self.no_sequences):
                    self.get_pdb(sequences_id[i])
                    c_file.writelines(
                        sequences_id[i] + ' ' + os.path.abspath('strike') + '/' + sequences_id[i] + '.pdb ' + chains[
                            i] + '\n')
                    a_file.writelines('>' + sequences_id[i] + '\n' + self.aligned_sequences[i] + '\n')

        return StrikeEx(os.path.abspath(self.exe_path)).run(
            parameters={'-c': self.out_connection_path, '-a': self.out_alignment_path})

    @staticmethod
    def get_pdb(pdb_id: str) -> str:
        base_url = 'https://files.rcsb.org/download/'
        base_url += pdb_id + '.pdb'

        out_pdb_path = 'strike/' + pdb_id + '.pdb'

        if not Path(out_pdb_path).is_file():
            LOGGER.debug('Downloading .pdb file...')

            try:
                req = urllib.request.Request(base_url)
                f = urllib.request.urlopen(req)
            except HTTPError:
                raise Exception("PDB not found in rcsb")

            with open(out_pdb_path, 'w+') as pdb_file:
                pdb_file.write(f.read().decode('unicode_escape'))

        return out_pdb_path

    @staticmethod
    def is_minimization() -> bool:
        return False
