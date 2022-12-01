"""
sequence.py
This file is part of the "Recherche de structure dans les protéines" project.
INFO-F-302 - Informatique fondamentale

Tue 29 Nov 2022

Authors:
    - Anton Romanova
    - Ismael Secundar

This file contains the protein sequence class.
"""
import functools

from src import validators


class Sequence:
    def __init__(self, sequence):
        validators.validate_sequence(sequence)

        self.sequence = sequence

        self._n = len(sequence)

    @property
    def n(self):
        return self._n

    @functools.cached_property
    def indices_of_1s(self):
        return tuple(i for i, v in enumerate(self.sequence) if v == "1")

    @property
    def number_of_1s(self):
        return len(self.indices_of_1s)

    @property
    def number_of_0s(self):
        return self.n - self.number_of_1s

    @property
    def is_all_ones(self):
        return self.number_of_1s == self.n

    @property
    def is_all_zeros(self):
        return self.number_of_0s == self.n

    @functools.lru_cache(maxsize=None)
    def get_flat_sequence_score(self):
        """
        Counts the number of adjacent 1s in the sequence.

        >>> Sequence("10111011011101").get_flat_sequence_score()
        5
        >>> Sequence("111").get_flat_sequence_score()
        2
        """
        return len([(i, j) for i, j in zip(self.sequence, self.sequence[1:]) if i == j == "1"])

    def __repr__(self):
        return f"Sequence({self.sequence!r})"

    def __str__(self):
        return self.sequence

    def __len__(self):
        return self.n

    def __getitem__(self, item):
        return self.sequence[item]

    def __iter__(self):
        return iter(self.sequence)

    def __eq__(self, other):
        return self.sequence == other.sequence

    def __hash__(self):
        return hash(self.sequence)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
