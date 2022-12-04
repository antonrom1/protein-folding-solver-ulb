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

    >>> validate_sequence('0')
    >>> validate_sequence('1')
    >>> validate_sequence('01')
    >>> validate_sequence('10')
    >>> validate_sequence('0110')
    >>> validate_sequence('1010')
    >>> validate_sequence('012')
    Traceback (most recent call last):
    ...
    ValidationError: The sequence must be non-empty and only contain '1' and '0'
    >>> validate_sequence('')
    Traceback (most recent call last):
    ...
    ValidationError: The sequence must be non-empty and only contain '1' and '0'
    >>> validate_sequence(1)
    Traceback (most recent call last):
    ...
    ValidationError: The sequence must be a string
    """
    if not isinstance(sequence, str):
        raise ValidationError("The sequence must be a string")

    if not re.match(r"^[01]+$", sequence):
        raise ValidationError("The sequence must be non-empty and only contain '1' and '0'")


def validate_bound(bound: int):
    """
    Validates the protein minimum score bound.

    :param bound: the protein minimum score bound (integer)

    >>> validate_bound(0)
    >>> validate_bound(1)
    >>> validate_bound('2')
    Traceback (most recent call last):
    ...
    ValidationError: The bound must be an integer
    >>> validate_bound(-1)
    Traceback (most recent call last):
    ...
    ValidationError: The bound must be positive
    """
    if not isinstance(bound, int):
        raise ValidationError("The bound must be an integer")

    if bound < 0:
        raise ValidationError("The bound must be positive")


if __name__ == "__main__":
    import doctest

    doctest.testmod()