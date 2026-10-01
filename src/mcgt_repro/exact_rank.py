"""Exact rank and nullspace over the rational numbers."""

from fractions import Fraction


def _frac(value):
    return value if isinstance(value, Fraction) else Fraction(value)


def rank_rows(rows):
    """Exact rank of rational rows, by Gauss--Jordan elimination."""
    if not rows:
        return 0
    matrix = [[_frac(entry) for entry in row] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = None
        for row in range(pivot_row, row_count):
            if matrix[row][column] != 0:
                pivot = row
                break
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pivot_value = matrix[pivot_row][column]
        matrix[pivot_row] = [entry / pivot_value for entry in matrix[pivot_row]]
        for row in range(row_count):
            if row != pivot_row and matrix[row][column] != 0:
                factor = matrix[row][column]
                matrix[row] = [
                    entry - factor * pivot_entry
                    for entry, pivot_entry in zip(matrix[row], matrix[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def nullspace(rows, column_count):
    """Exact basis of the rational kernel, as a list of column vectors."""
    if not rows:
        return [
            [Fraction(1) if i == j else Fraction(0) for i in range(column_count)]
            for j in range(column_count)
        ]
    matrix = [[_frac(entry) for entry in row] for row in rows]
    row_count = len(matrix)
    pivot_row = 0
    pivots = []
    for column in range(column_count):
        pivot = None
        for row in range(pivot_row, row_count):
            if matrix[row][column] != 0:
                pivot = row
                break
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pivot_value = matrix[pivot_row][column]
        matrix[pivot_row] = [entry / pivot_value for entry in matrix[pivot_row]]
        for row in range(row_count):
            if row != pivot_row and matrix[row][column] != 0:
                factor = matrix[row][column]
                matrix[row] = [
                    entry - factor * pivot_entry
                    for entry, pivot_entry in zip(matrix[row], matrix[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break
    pivot_set = set(pivots)
    basis = []
    for free_column in (c for c in range(column_count) if c not in pivot_set):
        vector = [Fraction(0)] * column_count
        vector[free_column] = Fraction(1)
        for index, pivot_column in enumerate(pivots):
            vector[pivot_column] = -matrix[index][free_column]
        basis.append(vector)
    return basis


def rank_columns(columns):
    if not columns:
        return 0
    return rank_rows([list(row) for row in zip(*columns)])


def apply_matrix(rows, columns):
    """Columns of the product of a row-stored matrix and column vectors."""
    images = []
    for column in columns:
        image = []
        for row in rows:
            total = Fraction(0)
            for left, right in zip(row, column):
                if left and right:
                    total += left * right
            image.append(total)
        images.append(image)
    return images


def is_zero_columns(columns):
    return all(entry == 0 for column in columns for entry in column)
