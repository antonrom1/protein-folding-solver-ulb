"""
folding.py
This file is part of the "Recherche de structure dans les protéines" project.
INFO-F-302 - Informatique fondamentale

Tue 29 Nov 2022

Authors:
    - Anton Romanova
    - Ismael Secundar

This file contains the protein folding solver.
"""
from __future__ import annotations

import itertools
import logging
import math
from dataclasses import dataclass
from textwrap import dedent
from typing import Optional, List

from src import sequence
from src.encoder import Encoder

logger = logging.getLogger(__name__)


class FoldedProtein:
    """
    FoldedProtein is a class that represents a folded protein.
    its main function is to be able to print itself in a nice way.
    """

    def __init__(self, folding: List[List[bool | None]]):
        self._folding = folding
        self._n = len(folding)
        assert all(len(row) == self._n for row in folding), "The folding must be a square matrix"

    def __repr__(self):
        return f"FoldedProtein({self._folding})"

    def __str__(self):
        output = ""

        for i, row in enumerate(self._folding):
            for j, v in enumerate(row):
                output += str(int(v)) if v is not None else " "
                output += ("—" if v and row[j + 1] else " ") if j + 1 < self._n else ""

            if i < self._n - 1:
                output += "\n"
                for j, v in enumerate(row):
                    output += ("|" if v and self._folding[i + 1][j] else " ") + (" " if j + 1 < self._n else "")
                output += "\n"

        # if we have empty lines in the beginning or end, remove them

        # reduce the empty spaces
        output = dedent(output)
        output = '\n'.join(map(str.rstrip, output.split('\n')))
        output = output.strip('\n')

        return output

    @classmethod
    def from_straight_sequence(cls, sequence: str | sequence.Sequence) -> FoldedProtein:
        """
        Create a FoldedProtein from a straight sequence.
        """
        first_row = [bool(int(v)) for v in sequence]
        remaining_rows = [[None] * len(first_row) for _ in range(len(first_row) - 1)]
        return cls([first_row] + remaining_rows)

    @classmethod
    def from_all_ones(cls, n: int) -> FoldedProtein:
        """
        Create a FoldedProtein from a sequence of all ones.
        """
        return cls.from_straight_sequence("1" * n)

    @classmethod
    def from_all_zeros(cls, n: int) -> FoldedProtein:
        """
        Create a FoldedProtein from a sequence of all zeros.
        """
        return cls.from_straight_sequence("0" * n)

    @classmethod
    def build_best_solution_for_all_ones(cls, n: int) -> FoldedProtein:
        """
        Build the best solution for a sequence of all ones.

        The solution is built according to the following rules specified here:
        https://www.irit.fr/~Jean-Michel.Hautiere/teaching/2021-2022/InfoFond/Projet/Projet.html#solution-optimale-pour-une-s
        https://oeis.org/A000041
        """
        matrix_size = math.ceil(n ** 0.5) + 1
        matrix = [[None] * matrix_size for _ in range(matrix_size)]

        # first we build a square with a side length of sqrt(n)
        biggest_square_that_fits = int(math.sqrt(n))
        for i in range(biggest_square_that_fits):
            for j in range(biggest_square_that_fits):
                matrix[i][j] = '1'

        remaining_ones = n - biggest_square_that_fits ** 2

        assert remaining_ones <= 2 * biggest_square_that_fits

        # then we add the remaining ones
        # first we add min(remaining_ones, biggest_square_that_fits) ones at the bottom

        for i in range(min(remaining_ones, biggest_square_that_fits)):
            matrix[biggest_square_that_fits][i] = '1'
            remaining_ones -= 1
            if remaining_ones == 0:
                break

        # then we add the remaining ones at the right
        for i in range(remaining_ones):
            matrix[i][biggest_square_that_fits] = '1'
            remaining_ones -= 1
            if remaining_ones == 0:
                break

        assert remaining_ones == 0

        return cls(matrix)


@dataclass
class FoldingSolution:
    """
    FoldingSolution is a class that represents a solution to the protein folding problem.
    """
    sequence: sequence.Sequence
    bound: int
    score: Optional[int]
    solution: Optional[FoldedProtein]

    @property
    def is_sat(self):
        return self.solution is not None

    @property
    def is_unsat(self):
        return self.solution is None


def solve_for_n(args):
    n, seq, bound = args
    logger.debug(f"Trying to solve for n={n}")
    m = seq.n - n + 1
    encoder = Encoder(seq, n, m, bound)
    return encoder.solve()


class FoldingSolver:
    """
    FoldingSolver is a class that solves the protein folding problem.
    In addition to the encoding the problem with Encoder,
    it also contains the logic to solve trivial cases in a more efficient way.
    """
    def __init__(self, seq: sequence.Sequence):
        self._sequence = seq
        self._max_possible_bound = self.sequence.max_score_bound

    @property
    def sequence(self):
        return self._sequence

    @property
    def max_possible_score(self):
        return self._max_possible_bound

    def solve(self, bound) -> FoldingSolution:
        if (res := self._test_for_trivial_cases(bound)) is not None:
            return res

        matrix_size = int(math.ceil(self.sequence.n ** (2 / 3)))
        # 1 + self.sequence.n // 4 if self.sequence.n >= 12 else self.sequence.n

        encoder = Encoder(self.sequence, matrix_size, matrix_size, bound)
        return encoder.solve()

    def _model_to_solution(self, model: List[int], bound: int, encoder: Encoder) -> FoldingSolution:
        """
        Convert a model to a solution.
        """
        if model is None:
            return FoldingSolution(self.sequence, bound, None, None)

        literals_set = set(model)
        protein_matrix = [[None] * self.sequence.n for _ in range(self.sequence.n)]

        for i, j in itertools.product(range(self.sequence.n - 2), repeat=2):
            for value_idx, value in enumerate(self.sequence):
                if encoder.get_literal(i, j, value_idx) in literals_set:
                    protein_matrix[i][j] = value if value is None else bool(int(value))
                    break

        n_bonds = sum(1 for bond in encoder.bond_literals if bond in literals_set)

        return FoldingSolution(self.sequence, bound, n_bonds, FoldedProtein(protein_matrix))

    def _test_for_trivial_cases(self, bound) -> Optional[FoldingSolution]:
        """
        Test for trivial cases where the sequence is already a solution.
        """
        if self.sequence.is_all_ones:
            if bound <= self.max_possible_score:
                logger.info("Trivial solution: Sequence is all ones, and bound is less than max possible score.")
                # assert False, "Bad solution"
                return FoldingSolution(self.sequence, bound, self.max_possible_score,
                                       FoldedProtein.build_best_solution_for_all_ones(self.sequence.n))
        elif self.sequence.is_all_zeros:
            if bound == 0:
                logger.info("Trivial solution: Sequence is all zeros, and bound is 0.")
                return FoldingSolution(self.sequence, bound, 0, FoldedProtein.from_all_zeros(self.sequence.n))

        assert (flat_sequence_score := self.sequence.get_flat_sequence_score()) <= self.max_possible_score

        if bound > self.max_possible_score:
            logger.info("Trivial solution: Bound is greater than max possible score.")
            return FoldingSolution(self.sequence, bound, None, None)

        if flat_sequence_score >= bound:
            logger.info("Trivial solution: Flat sequence score is greater than bound.")
            return FoldingSolution(self.sequence, bound, flat_sequence_score,
                                   FoldedProtein.from_straight_sequence(self.sequence.sequence))


if __name__ == "__main__":
    import doctest

    doctest.testmod()
