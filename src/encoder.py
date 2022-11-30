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

from pysat.card import CNF, IDPool, CardEnc

from src import sequence

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
    def __init__(self, seq: sequence.Sequence, bound):
        self.sequence = seq
        self.bound = bound

        self._cnf = CNF()
        self._idp = IDPool()

        self._sequence_indices = set(range(len(self.sequence)))
        self._matrix_variables = self._sequence_indices | {None}

        self._init_grid_constraints()
        self._cnf.extend(self._add_bound_constraint())

        logger.info(f"initialized encoder with {len(self._cnf.clauses)} clauses")

    @property
    def cnf(self):
        return self._cnf

    def get_literal(self, i, j, value_idx):
        return self._idp.id((i, j, value_idx))

    def _init_grid_constraints(self):
        self._cnf.extend(self._one_value_per_position())
        self._cnf.extend(self._unique_value_in_matrix())
        self._cnf.extend(self._sequence_values_adjacent())

    @clauses_logging_description("every matrix position has only one value")
    def _one_value_per_position(self):
        clauses = []
        for i, j in itertools.product(range(len(self.sequence)), repeat=2):
            possible_literals = [self._idp.id((i, j, v)) for v in self._matrix_variables]
            clauses.extend(CardEnc.equals(possible_literals, 1, vpool=self._idp))
        return clauses

    @clauses_logging_description("every value (except None) appears only once in the matrix")
    def _unique_value_in_matrix(self):
        clauses = []
        for v in self._sequence_indices:
            clauses.extend(
                CardEnc.equals(
                    [self._idp.id((i, j, v)) for i, j in itertools.product(range(len(self.sequence)), repeat=2)],
                    1,
                    vpool=self._idp
                )
            )
        return clauses

    @clauses_logging_description("every value is adjacent to its neighbours in the sequence")
    def _sequence_values_adjacent(self):
        clauses = []
        for i, j in itertools.product(range(len(self.sequence)), repeat=2):
            # TODO: this can be optimized by skipping 1/2 of the values because they are already checked
            for v in range(1, len(self.sequence) - 1):
                # [(i, j, v) => [(i - 1, j, v - 1) | (i + 1, j, v - 1) | (i, j - 1, v - 1) | (i, j + 1, v - 1)]]
                # a => (b | c | d | e) <=> !a | b | c | d | e

                v1p1_clause = [-self._idp.id((i, j, v))]
                v1m1_clause = [-self._idp.id((i, j, v))]

                for i2, j2 in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
                    if 0 <= i2 < len(self.sequence) and 0 <= j2 < len(self.sequence):
                        v1p1_clause.append(self._idp.id((i2, j2, v - 1)))
                        v1m1_clause.append(self._idp.id((i2, j2, v + 1)))
                clauses.append(v1m1_clause)
                clauses.append(v1p1_clause)

        return clauses

    @clauses_logging_description("there are at least 'bound' adjacent 1s")
    def _add_bound_constraint(self):
        clauses = []
        clauses_for_at_least = []
        indices_of_1s = [i for i, v in enumerate(self.sequence) if v == '1']
        for i, j in itertools.product(range(len(self.sequence)), repeat=2):

            if (i + j) % 2 == 0:
                # this is a little optimization to avoid checking the same clauses twice
                continue

            for v1, v2 in itertools.permutations(indices_of_1s, 2):
                for i2, j2 in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
                    if not (0 <= i2 < len(self.sequence) and 0 <= j2 < len(self.sequence)):
                        continue
                    a = self._idp.id((i, j, v1))
                    b = self._idp.id((i2, j2, v2))
                    c = self._idp.id(((i, i2), (j, j2), (v1, v2)))

                    # (a & b) <=> c
                    # cnf: (¬c ∨ a) ∧ (¬c ∨ b) ∧ (¬a ∨ ¬b ∨ c)
                    clauses.extend([[-c, a], [-c, b], [a, b, -c]])
                    clauses_for_at_least.append(c)

        clauses.extend(CardEnc.atleast(clauses_for_at_least, self.bound, vpool=self._idp))
        return clauses
