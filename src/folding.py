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
import re
from dataclasses import dataclass
from typing import Optional, List

from pysat.solvers import Glucose4

from src import sequence, utils
from src.encoder import Encoder


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

        output = output.splitlines()

        # remove empty lines in the beginning
        while output and re.match(r'^\s*$', output[0]):
            output.pop(0)

        # remove empty lines in the end
        while output and re.match(r'^\s*$', output[-1]):
            output.pop()

        output = "\n".join(output)

        # remove empty rows on the left and right
        # find min amount of spaces on the left and right
        min_spaces_left = min(len(row) - len(row.lstrip()) for row in output.splitlines())
        min_spaces_right = min(len(row) - len(row.rstrip()) for row in output.splitlines())

        # remove the spaces
        output = "\n".join(
            row[min_spaces_left:(-min_spaces_right) if min_spaces_right else None] for row in output.splitlines())

        return output


    @classmethod
    def from_straight_sequence(cls, sequence: str) -> FoldedProtein:
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


@dataclass
class FoldingSolution:
    sequence: sequence.Sequence
    bound: int
    solution: Optional[FoldedProtein]

    @property
    def is_sat(self):
        return self.solution is not None

    @property
    def is_unsat(self):
        return self.solution is None


class FoldingSolver:
    def __init__(self, seq: sequence.Sequence):
        self._sequence = seq
        self._max_possible_score = utils.get_max_score_bound_for_length(self.sequence.n)

    @property
    def sequence(self):
        return self._sequence

    @property
    def max_possible_score(self):
        return self._max_possible_score

    def solve(self, bound) -> FoldingSolution:
        if (res := self._test_for_trivial_cases(bound)) is not None:
            return res

        encoder = Encoder(self.sequence, bound)
        cnf = encoder.cnf

        solver = Glucose4(use_timer=True)
        solver.append_formula(cnf)
        solver.solve()

        model = solver.get_model()

        return self._model_to_solution(model, bound, encoder)

    def _model_to_solution(self, model: List[int], bound: int, encoder: Encoder) -> FoldingSolution:
        """
        Convert a model to a solution.
        """
        if model is None:
            return FoldingSolution(self.sequence, bound, None)

        literals_set = set(model)
        protein_matrix = [[None] * self.sequence.n for _ in range(self.sequence.n)]

        for i, j in itertools.product(range(self.sequence.n), repeat=2):
            for value_idx, value in enumerate(self.sequence):
                if encoder.get_literal(i, j, value_idx) in literals_set:
                    protein_matrix[i][j] = value if value is None else bool(int(value))
                    break

        return FoldingSolution(self.sequence, bound, FoldedProtein(protein_matrix))



    def _test_for_trivial_cases(self, bound) -> Optional[FoldingSolution]:
        """
        Test for trivial cases where the sequence is already a solution.

        >>> FoldingSolver(sequence.Sequence("1"))._test_for_trivial_cases(0)
        FoldingSolution(sequence=Sequence('1'), bound=0, solution=FoldedProtein([[True]]))

        >>> FoldingSolver(sequence.Sequence("1"))._test_for_trivial_cases(1)
        FoldingSolution(sequence=Sequence('1'), bound=1, solution=None)

        >>> FoldingSolver(sequence.Sequence("0"))._test_for_trivial_cases(0)
        FoldingSolution(sequence=Sequence('0'), bound=0, solution=FoldedProtein([[False]]))

        >>> FoldingSolver(sequence.Sequence('1' * 4))._test_for_trivial_cases(4)
        FoldingSolution(sequence=Sequence('1111'), bound=4, solution=FoldedProtein([[True, True, True, True], [None, None, None, None], [None, None, None, None], [None, None, None, None]]))

        """
        if self.sequence.is_all_ones:
            if bound <= self.max_possible_score:
                return FoldingSolution(self.sequence, bound, FoldedProtein.from_all_ones(self.sequence.n))
        elif self.sequence.is_all_zeros:
            if bound == 0:
                return FoldingSolution(self.sequence, bound, FoldedProtein.from_all_zeros(self.sequence.n))
            else:
                return FoldingSolution(self.sequence, bound, None)

        if bound >= self.max_possible_score:
            # only full ones sequences can have a score of max_possible_score
            return FoldingSolution(self.sequence, bound, None)

        if self.sequence.get_flat_sequence_score() >= bound:
            return FoldingSolution(self.sequence, bound, FoldedProtein.from_straight_sequence(self.sequence.sequence))


if __name__ == "__main__":
    import doctest

    doctest.testmod()

    s = sequence.Sequence('1110101101')
    solver = FoldingSolver(s)
    print(solver.solve(7).solution)
