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
import math

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

    @property
    def max_score_bound(self):
        """
        Get the maximum score bound for a sequence of length n.
        The maximum score bound can only be reached if the sequence is all 1s and arranged as such:

        If one constructs a square (square 1) and then draws another square of identical size beside it (square 2),
        the squares share 1 edge. If one then places an identical square above square 2
        (instead of continuing in a straight path), there are now 2 shared edges.
        Continuing this pattern in an outward spiral, one finds that the number of shared edges is 4, 5, 7,

        For more information, see https://oeis.org/A123663.

        :return: the maximum score bound

        >>> [Sequence('1' * i).max_score_bound for i in range(1, 26)]
        [0, 1, 2, 4, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 24, 25, 27, 29, 31, 32, 34, 36, 38, 40]
        """

        # implementation by Chai Wah Wu, Jul 28 2022, https://oeis.org/A123663
        return (m := self.number_of_1s << 1) - 1 - math.isqrt((m << 1) - 1) if self.number_of_1s > 0 else 0

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
