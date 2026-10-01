"""Real bases of Herm(r), u(r), and Mat(r, C), with rational coordinates."""

from fractions import Fraction


class QC(object):
    __slots__ = ("re", "im")

    def __init__(self, real=0, imag=0):
        self.re = real if isinstance(real, Fraction) else Fraction(real)
        self.im = imag if isinstance(imag, Fraction) else Fraction(imag)

    def __add__(self, other):
        return QC(self.re + other.re, self.im + other.im)

    def __sub__(self, other):
        return QC(self.re - other.re, self.im - other.im)

    def __neg__(self):
        return QC(-self.re, -self.im)

    def __mul__(self, other):
        if isinstance(other, QC):
            return QC(
                self.re * other.re - self.im * other.im,
                self.re * other.im + self.im * other.re,
            )
        scale = other if isinstance(other, Fraction) else Fraction(other)
        return QC(self.re * scale, self.im * scale)

    def conj(self):
        return QC(self.re, -self.im)

    def __eq__(self, other):
        return self.re == other.re and self.im == other.im


ZERO = QC(0, 0)
ONE = QC(1, 0)
I_UNIT = QC(0, 1)


def zeros(size):
    return tuple(tuple(ZERO for _ in range(size)) for _ in range(size))


def add(left, right):
    return tuple(
        tuple(a + b for a, b in zip(row_left, row_right))
        for row_left, row_right in zip(left, right)
    )


def sub(left, right):
    return tuple(
        tuple(a - b for a, b in zip(row_left, row_right))
        for row_left, row_right in zip(left, right)
    )


def scale(matrix, factor):
    return tuple(tuple(entry * factor for entry in row) for row in matrix)


def adjoint(matrix):
    size = len(matrix)
    return tuple(
        tuple(matrix[i][j].conj() for i in range(size)) for j in range(size)
    )


def herm_part(matrix):
    return scale(add(matrix, adjoint(matrix)), Fraction(1, 2))


def skew_part(matrix):
    return scale(sub(matrix, adjoint(matrix)), Fraction(1, 2))


def _unit(size, row, column, value):
    return tuple(
        tuple(value if (i, j) == (row, column) else ZERO for j in range(size))
        for i in range(size)
    )


def herm_basis(sector):
    """Real basis of Hermitian matrices. Its length is sector squared."""
    basis = []
    for index in range(sector):
        basis.append(_unit(sector, index, index, ONE))
    for row in range(sector):
        for column in range(row + 1, sector):
            basis.append(add(_unit(sector, row, column, ONE), _unit(sector, column, row, ONE)))
            basis.append(
                add(_unit(sector, row, column, I_UNIT), _unit(sector, column, row, QC(0, -1)))
            )
    return basis


def skew_basis(sector):
    """Real basis of skew-Hermitian matrices. Its length is sector squared."""
    basis = []
    for index in range(sector):
        basis.append(_unit(sector, index, index, I_UNIT))
    for row in range(sector):
        for column in range(row + 1, sector):
            basis.append(
                add(_unit(sector, row, column, ONE), _unit(sector, column, row, QC(-1, 0)))
            )
            basis.append(
                add(_unit(sector, row, column, I_UNIT), _unit(sector, column, row, I_UNIT))
            )
    return basis


def gl_basis(sector):
    """Real basis of Mat(sector, C), Hermitian part followed by skew-Hermitian part."""
    return herm_basis(sector) + skew_basis(sector)


def herm_coords(matrix):
    size = len(matrix)
    coordinates = [matrix[index][index].re for index in range(size)]
    for row in range(size):
        for column in range(row + 1, size):
            coordinates.append(matrix[row][column].re)
            coordinates.append(matrix[row][column].im)
    return coordinates


def skew_coords(matrix):
    size = len(matrix)
    coordinates = [matrix[index][index].im for index in range(size)]
    for row in range(size):
        for column in range(row + 1, size):
            coordinates.append(matrix[row][column].re)
            coordinates.append(matrix[row][column].im)
    return coordinates


def gl_coords(matrix):
    return herm_coords(herm_part(matrix)) + skew_coords(skew_part(matrix))
