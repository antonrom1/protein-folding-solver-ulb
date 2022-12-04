"""
validators.py
This file is part of the "Recherche de structure dans les protéines" project.
INFO-F-302 - Informatique fondamentale

Tue 29 Nov 2022

Authors:
    - Anton Romanova
    - Ismael Secundar

This file contains sequence validators.
"""

import re


class ValidationError(Exception):
    pass


def validate_sequence(sequence: str):
    """
    Validate the protein sequence.

    :param sequence: the protein sequence (string)
    """
    if not isinstance(sequence, str):
        raise ValidationError("The sequence must be a string")

    if not re.match(r"^[01]+$", sequence):
        raise ValidationError("The sequence must be non-empty and only contain '1' and '0'")


def validate_bound(bound: int):
    """
    Validates the protein minimum score bound.

    :param bound: the protein minimum score bound (integer)
    """
    if not isinstance(bound, int):
        raise ValidationError("The bound must be an integer")

    if bound < 0:
        raise ValidationError("The bound must be positive")
