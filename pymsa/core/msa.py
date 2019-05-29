from typing import List


class MSA:

    def __init__(self, sequences: list, ids: list = None, gap_character: str = '-'):
        self._sequences = sequences
        self._ids = ids
        self.gap_character = gap_character

    @property
    def sequences(self) -> List[str]:
        """
        :return: Aligned sequences.
        """
        return self._sequences

    @property
    def ids(self):
        """
        :return: Sequences identifiers.
        """
        return self._ids

    @property
    def number_of_sequences(self) -> int:
        """
        :return: Number of sequences within the alignment.
        """
        return len(self._sequences)

    @property
    def is_valid(self) -> bool:
        return all(len(seq) == len(self.sequences[0]) for seq in self._sequences[:1]) and len(self._sequences) >= 2

    def __len__(self) -> int:
        """
        :return: Total length of the alignment.
        """
        return len(self.sequences[0])
