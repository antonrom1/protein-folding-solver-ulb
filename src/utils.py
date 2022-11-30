"""
utils.py
This file is part of the "Recherche de structure dans les protéines" project.
INFO-F-302 - Informatique fondamentale

Tue 29 Nov 2022

Authors:
    - Anton Romanova
    - Ismael Secundar

This file contains utility functions.
"""

import math


def get_max_score_bound_for_length(n):
    """
    Get the maximum score bound for a sequence of length n.
    The maximum score bound can only be reached if the sequence is all 1s and arranged as such:

    If one constructs a square (square 1) and then draws another square of identical size beside it (square 2),
    the squares share 1 edge. If one then places an identical square above square 2
    (instead of continuing in a straight path), there are now 2 shared edges.
    Continuing this pattern in an outward spiral, one finds that the number of shared edges is 4, 5, 7,

    For more information, see https://oeis.org/A000041.

    :param n: the length of the sequence
    :return: the maximum score bound

    >>> list(map(get_max_score_bound_for_length, range(1, 26)))
    [0, 1, 2, 4, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 24, 25, 27, 29, 31, 32, 34, 36, 38, 40]
    """
    return (m := n << 1) - 1 - math.isqrt((m << 1) - 1)


if __name__ == '__main__':
    import doctest

    doctest.testmod()
