from typing import *
import fractions as fr


class Matrix:
    """A matrix like in linear algebra
    _elements: the elements in a matrix denoted as a list of lists.
    _augmented: whether the last column is augmented or not.
    Follows Matrix[row_num][col_num]
    """
    _elements: List[List[fr.Fraction]]
    _augmented: bool

    def __init__(self, elements, augmented=False) -> None:
        """Precondition: elements represents a valid matrix"""
        for sublist in elements:
            for i in range(len(sublist)):
                sublist[i] = fr.Fraction(sublist[i]).limit_denominator()

        self._elements = elements
        self._augmented = augmented

    def __str__(self) -> str:
        if not self._augmented:
            s = ""
            for row in self._elements:
                s += '['
                for element in row:
                    s += "{} ".format(str(element))
                s = s.strip()
                s += ']\n'
            return s
        else:
            s = ""
            for row in self._elements:
                s += '['
                for element in row[:len(row)-1]:
                    s += "{} ".format(str(element))
                s += "| {}]\n".format(str(row[-1]))
            return s

    def row_reduce(self) -> None:
        """Mutates this matrix into its row reduced echelon form
        Augmented determines whether the last column of the matrix is considered
        the constant part of a system of linear equations. True for yes, False
        for no."""
        non_zero_col = self._find_first_non_zero_column()
        if non_zero_col == -1:
            return

        # Track leading ones for converting REF to RREF
        leading_ones = []

        # Convert to row echelon form
        row_offset = 0
        if not self._augmented:
            last_column = len(self._elements[0])
        else:
            last_column = len(self._elements[0]) - 1

        for col_num in range(non_zero_col, last_column):
            if not self.check_zero_column_below_row_num(col_num, row_offset):
                if self._elements[row_offset][col_num] == 0:
                    self.swap_row(row_offset,
                                  self._find_non_zero(row_offset, col_num))
                scalar = 1 / self._elements[row_offset][col_num]
                leading_ones.append((row_offset, col_num))
                self.scale_row(row_offset, scalar)
                for row_num in range(row_offset + 1, len(self._elements)):
                    self.add_row(row_num, row_offset,
                                 self._elements[row_num][col_num] * -1)
                row_offset += 1

        # Convert to row reduced echelon form
        # leading_one: [(rowPosition, colPosition),(rowPosition, colPosition)]
        leading_ones.reverse()
        for leading_one in leading_ones:
            for row_num in range(leading_one[0]-1, -1, -1):
                scalar = self._elements[row_num][leading_one[1]] * -1
                self.add_row(row_num, leading_one[0], scalar)

    def check_zero_column_below_row_num(self, col_num, row_num) -> bool:
        """Returns whether all elements in column col_num below and
        including row row_num are 0"""
        for row_number in range(row_num, len(self._elements)):
            if self._elements[row_number][col_num] != 0:
                return False
        return True

    def _find_non_zero(self, row_offset: int, col_num: int) -> int:
        """returns the position of the closest row below offset that has
        a non-zero entry in column col_num. Returns -1 if no such element
        exists"""
        for row_num in range(row_offset, len(self._elements)):
            if self._elements[row_num][col_num] != 0:
                return row_num
        return -1

    def _find_first_non_zero_column(self) -> int:
        """Returns the position of the first non-zero column or -1 is the
        matrix is the zero matrix"""
        for col_num in range(len(self._elements[0])):
            for row_num in range(len(self._elements)):
                if self._elements[row_num][col_num] != 0:
                    return col_num
        return -1

    def swap_row(self, row_num_one: int, row_num_two: int) -> None:
        self._elements[row_num_one], self._elements[row_num_two] \
            = self._elements[row_num_two], self._elements[row_num_one]

    def scale_row(self, row_num: int, scalar: float) -> None:
        for i in range(len(self._elements[row_num])):
            self._elements[row_num][i] = self._elements[row_num][i] * scalar

    def add_row(self, row_num_one: int, row_num_two: int, scalar: float) -> \
            None:
        """Adds all the elements in row_num_two times scalar to the respective
        elements in row_num_one"""
        for i in range(len(self._elements[row_num_one])):
            self._elements[row_num_one][i] += self._elements[row_num_two][i] \
                                              * scalar
