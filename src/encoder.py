"""
encoder.py
This file is part of the "Recherche de structure dans les protéines" project.
INFO-F-302 - Informatique fondamentale

Tue 29 Nov 2022

Authors:
    - Anton Romanova
    - Ismael Secundar

This file contains the SAT encoder class.
"""
import itertools
import logging

from pysat.card import CNF, IDPool, CardEnc, EncType
from pysat.solvers import MapleCM, Glucose4

from src import sequence, folding

logger = logging.getLogger(__name__)


def clauses_logging_description(condition_description):
    def decorator(func):
        def wrapper(*args, **kwargs):
            logger.info(f'adding clauses to check if {condition_description}')
            result = func(*args, **kwargs)
            logger.debug(f"added {len(result)} clauses to check whether {condition_description}")
            return result

        return wrapper

    return decorator


class Encoder:
    def __init__(self, seq: sequence.Sequence, n, m, bound):
        self.sequence = seq
        self.bound = bound

        # we do not need to create a matrix of size n x n because
        # if the flat sequence already has the necessary score, the trivial cases checks will be enough
        self._n = n
        self._m = m

        # assert self._n + self._m - 1 <= self.sequence.n, "The matrix is too big for the sequence"
        # assert self._n * self._m >= self.sequence.n, "The matrix is too small for the sequence"

        self._bond_literals = []

        self._cnf = CNF()
        self._idp = IDPool()

        self._sequence_indices = set(range(len(self.sequence)))
        self._matrix_variables = self._sequence_indices | {None}
        self._matrix_positions = set(itertools.product(range(self.n), range(self.m)))

        self._init_grid_constraints()
        self._cnf.extend(self._add_at_least_bound_bonds_constraint())

        logger.info(f"initialized encoder with {len(self._cnf.clauses)} clauses")

    def solve(self):
        solver = Glucose4(bootstrap_with=self.cnf.clauses)
        solver.solve()
        return self._model_to_solution(solver.get_model())

    def _model_to_solution(self, model) -> 'folding.FoldingSolution':
        """
        Convert a model to a solution.
        """
        if model is None:
            return folding.FoldingSolution(self.sequence, self.bound, None, None)

        literals_set = set(model)
        protein_matrix = [[None] * self.sequence.n for _ in range(self.sequence.n)]

        for i, j in itertools.product(range(self.sequence.n - 2), repeat=2):
            for value_idx, value in enumerate(self.sequence):
                if self.get_literal(i, j, value_idx) in literals_set:
                    protein_matrix[i][j] = value if value is None else bool(int(value))
                    break

        n_bonds = sum(1 for bond in self.bond_literals if bond in literals_set)

        return folding.FoldingSolution(self.sequence, self.bound, n_bonds, folding.FoldedProtein(protein_matrix))

    @property
    def n(self):
        return self._n

    @property
    def m(self):
        return self._m

    @property
    def cnf(self):
        return self._cnf

    @property
    def bond_literals(self):
        return tuple(self._bond_literals)

    def get_literal(self, i, j, value_idx):
        return self._idp.id((i, j, value_idx))

    def _init_grid_constraints(self):
        self._cnf.extend(self._one_value_per_position())
        self._cnf.extend(self._unique_value_in_matrix())
        self._cnf.extend(self._sequence_values_adjacent())

    @clauses_logging_description("every matrix position has only one value")
    def _one_value_per_position(self):
        clauses = []
        for i, j in self._matrix_positions:
            clauses.append([self._idp.id((i, j, v)) for v in self._matrix_variables])

        # at most one value per position
        for i, j in self._matrix_positions:
            for v1, v2 in itertools.combinations(self._matrix_variables, 2):
                clauses.append([-self._idp.id((i, j, v1)), -self._idp.id((i, j, v2))])
        return clauses

    @clauses_logging_description("every value (except None) appears only once in the matrix")
    def _unique_value_in_matrix(self):
        clauses = []
        for v in self._sequence_indices:
            clauses.extend(
                CardEnc.equals(
                    [self._idp.id((i, j, v)) for i, j in self._matrix_positions],
                    1,
                    vpool=self._idp
                )
            )
        return clauses

    @clauses_logging_description("every value is adjacent to its neighbours in the sequence")
    def _sequence_values_adjacent(self):
        clauses = []
        for i, j in self._matrix_positions:
            for v in range(1, len(self.sequence) - 1):
                # [(i, j, v) => [(i - 1, j, v - 1) | (i + 1, j, v - 1) | (i, j - 1, v - 1) | (i, j + 1, v - 1)]]
                # a => (b | c | d | e) <=> !a | b | c | d | e

                v1p1_clause = [-self._idp.id((i, j, v))]
                v1m1_clause = [-self._idp.id((i, j, v))]

                for i2, j2 in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
                    if 0 <= i2 < self._n and 0 <= j2 < self._m:
                        v1p1_clause.append(self._idp.id((i2, j2, v - 1)))
                        v1m1_clause.append(self._idp.id((i2, j2, v + 1)))
                clauses.append(v1m1_clause)
                clauses.append(v1p1_clause)

        return clauses

    @clauses_logging_description("there are at least 'bound' bonds")
    def _add_at_least_bound_bonds_constraint(self):
        clauses = []
        for i, j in self._matrix_positions:

            if (i + j) % 2 == 0:
                # avoid checking the same clauses twice
                continue

            adjacent_positions = filter(
                lambda pos: 0 <= pos[0] < self._n and 0 <= pos[1] < self._m,
                ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1))
            )

            for i2, j2 in adjacent_positions:
                for v1, v2 in itertools.permutations(self.sequence.indices_of_1s, 2):
                    cell = self._idp.id((i, j, v1))
                    adjacent_cell = self._idp.id((i2, j2, v2))
                    bond = self._idp.id(((i, i2), (j, j2), (v1, v2)))
                    # (cell & adjacent_cell) <=> bond
                    clauses.extend([[-bond, cell], [-bond, adjacent_cell]])
                    self._bond_literals.append(bond)
        clauses.extend(CardEnc.atleast(self._bond_literals, self.bound, vpool=self._idp))
        return clauses
