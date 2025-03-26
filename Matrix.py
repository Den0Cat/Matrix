from random import randint
from fractions import Fraction
from copy import deepcopy
import time    

def runtime(function):
    start = time.time()
    function()
    end = time.time()
    print(f"PROGRAM ENDS IN {end-start} sec.")


class Matrix():
    
    # get column space of the matrix
    def T(self):
        matrix: list[list[Fraction]] = deepcopy(self.matrix)
        col_space = []
        for j in range(len(matrix[0]) if isinstance(matrix[0], list) else len(matrix)):
            col_space.append([])
            for i in range(len(matrix) if isinstance(matrix[0], list) else 1):
                col_space[j].append(matrix[i][j] if isinstance(matrix[0], list) else matrix[j])
        return Matrix(col_space, f"({self.name})^T") 
    
    # convert matrix to Reduced Row Echelon Form (RREF)
    def rref(self):
        matrix: list[list[Fraction]] = deepcopy(self.matrix)
        pivots = dict()
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                if i not in pivots and matrix[i][j]!=0:
                    pivots[i] = matrix[i][j], i, j
                if i in pivots:
                    pivot = pivots[i][0]
                    for piv_i in range(len(matrix)):
                        if piv_i!=i:
                            k = (-matrix[piv_i][j]) / pivot
                            for piv_j in range(j, len(matrix[0])):
                                matrix[piv_i][piv_j] += matrix[i][piv_j]*k
                    break
        for i in pivots:
            for j in range(len(matrix[0])):
                matrix[i][j] /= pivots[i][0]
        return Matrix(matrix, f"({self.name})_rref"), pivots

    # get the null space from matrix
    def nullspace(self, printable: bool = False):
        RREF, pivots = deepcopy(self.rref())
        x = [f"x_{i+1}" for i in range(len(RREF.matrix[0]))]
        if printable:
            print("To find the null space, we need to solve the equation Ax=0 or Rx=0\n")
            print(f"x = {x}^T")
            print("R (Reduced Row Echelon Form):")
            RREF.print_matrix()
            print()
        pivot_vars = [pivots[i][2] for i in pivots]
        free_vars = [j for j in range(len(x)) if j not in pivot_vars]
        pivot_i = [pivots[i][1] for i in pivots]
        if printable:
            print("Pivot variables:", *[x[i] for i in pivot_vars])
            print("Free variables:", *[x[i] for i in range(len(x)) if i not in pivot_vars])
            print()
        if not free_vars:
            if printable:
                NullSpace = Matrix([[0] * len(pivot_vars)], "NULL SPACE")
                print(NullSpace.name)
                NullSpace.print_matrix(col_space=True)
            return Matrix([[0] * len(pivot_vars)], f"({self.name})x=0")
        null_basis = [[1 if j==i else 0 for j in range(len(x))] for i in free_vars]
        cur_pivot_i = 0
        for i in range(len(pivot_i)):
            if printable:
                print(f"{x[pivot_vars[cur_pivot_i]]} = ", end='')
            first_sign = True
            for j in range(len(RREF.matrix[0])):
                if j!=pivot_vars[cur_pivot_i] and RREF.matrix[pivot_i[i]][j]!=0:
                    if printable:
                        print(f"{'+' if RREF.matrix[pivot_i[i]][j]<0 and not first_sign else '-' if RREF.matrix[pivot_i[i]][j]>=0 else ''}{' ' if not first_sign else ''}{abs(-RREF.matrix[pivot_i[i]][j])}*{x[j]}", end=' ')
                    first_sign = False
                    null_basis[free_vars.index(j)][pivot_vars[cur_pivot_i]] = -RREF.matrix[pivot_i[i]][j]
            cur_pivot_i += 1
            if printable:
                print()
        if printable:
            print()
            NullSpace = Matrix(null_basis, "NULL SPACE")
            print(NullSpace.name)
            NullSpace.print_matrix(col_space=True)
        return Matrix(null_basis, f"({self.name})_nullspace")
    
    # get the left null space from matrix
    def left_nullspace(self, printable: bool = False):
        return self.T().nullspace(printable=printable)
    
    # is matrix NxN
    def isSquare(self) -> bool:
        matrix: list[list[Fraction]] = deepcopy(self.matrix)
        for line in matrix:
            if len(matrix) != len(line):
                return False
        return True
    
    # is matrix NxM
    def isLegal(self) -> bool:
        matrix: list[list[Fraction]] = deepcopy(self.matrix)
        if not (type(self.matrix[0]) is list):
            return True
        for i in range(len(matrix)):
            if len(matrix[0]) != len(matrix[i]):
                return False
        return True
    
    # inverse of the matrix
    def inverse(self):
        matrix: list[list[Fraction]] = deepcopy(self.matrix)
        if not self.isSquare():
            raise IndexError("MATRIX SHOULD BE SQUARE TO INVERSE")
        I = [[1 if i==j else 0 for j in range(len(matrix))] for i in range(len(matrix))]
        pivots = 1
        for j in range(len(matrix)):
            cur = matrix[j][j]
            if cur == 0:
                raise ValueError("DETERMINANT SHOULD BE NON ZERO TO INVERSE")
            pivots *= cur
            for i in range(len(matrix)):
                if i!=j:
                    k = -matrix[i][j]/cur
                    for t in range(j, len(matrix)):
                        matrix[i][t] += matrix[j][t]*k
                    for t in range(len(I)):
                        I[i][t] += I[j][t]*k
        for i in range(len(I)):
            cur = matrix[i][i]
            matrix[i][i] = cur / cur
            for j in range(len(I)):
                I[i][j] /= cur

        return Matrix(I, f"({self.name})^-1")
    
    # determinant of the matrix
    def det(self) -> Fraction:
        matrix: list[list[Fraction]] = deepcopy(self.matrix)
        I = [[1 if i==j else 0 for j in range(len(matrix))] for i in range(len(matrix))]
        pivots = 1
        for j in range(len(matrix)):
            cur = matrix[j][j]
            if cur == 0: return 0
            pivots *= cur
            for i in range(len(matrix)):
                if i!=j:
                    k = -matrix[i][j]/cur
                    for t in range(j, len(matrix)):
                        matrix[i][t] += matrix[j][t]*k
                    for t in range(len(I)):
                        I[i][t] += I[j][t]*k

        return pivots
    
    # length of the given vector
    def vec_len(self) -> Fraction:
        if isinstance(self.matrix[0], list) and len(self.matrix[0]) > 1:
            raise ValueError("LENGTH IS ONLY FOR VECTORS, NOT MATRICES")
        vec = deepcopy(self)
        if isinstance(self.matrix[0], list) and len(self.matrix[0]) == 1: vec = vec.T()
        length = Fraction(0)
        vec = vec.matrix
        for i in range(len(vec)):
            length += (vec[i])**2
        return length**0.5

    # orthonormal basis of the given matrix
    def orthonormal_basis(self):
        Q = deepcopy(self.T().matrix)
        for i in range(1, len(Q)):
            el = Matrix(Q[i], "element")
            fixed_el = deepcopy(el)
            for j in range(i):
                el_before = Matrix(Q[j], "before")
                el -= el_before * ( (fixed_el*el_before.T()).matrix[0] / (el_before*el_before.T()).matrix[0] )
            Q[i] = el.matrix
        for i in range(len(Q)):
            Q[i] = (Matrix(Q[i], "Q") / Matrix(Q[i], "Q").vec_len()).matrix
        Q = Matrix(Q, f"Q({self.name})").T()
        Q.name = f"Q({self.name})"
        return Q

    # A = QR; Q - orthonormal basis
    def Gram_Schimdt(self):
        Q = self.orthonormal_basis()
        R = Q.T()*self
        R.name = f"R({self.name})"
        return Q, R

    # text to array
    def parse_matrix(self, text: str) -> list[list[Fraction]]:
        if text.strip() == '':
            return ValueError("MATRIX FROM AN EMPTY STRING IS ILLEGAL")
        matrix_txt = text.split("\n")
        matrix: list = []
        for i in range(len(matrix_txt)):
            if matrix_txt[i].strip() != '':
                new_line = matrix_txt[i].split()
                matrix.append([Fraction(new_line[j]) for j in range(len(new_line))])
        return matrix
    
    # print matrix to the output stream
    def print_matrix(self, col_space: bool = False, print_name: bool = True):
        if print_name: print(self.name)
        matrix: list[list[Fraction]] = deepcopy(self.matrix) if not col_space else self.T().matrix
        print('', *[f"[{i}]" for i in range(len(matrix[0] if type(matrix[0]) is list else matrix))], sep='\t')
        if not isinstance(matrix[0], list): print("[0]", end='\t')
        for i in range(len(matrix)):
            if type(matrix[0]) is list:
                print(f"[{i}]", end='\t')
                for j in range(len(matrix[i])):
                    item = matrix[i][j]
                    if '/' in str(item):
                        a, b = str(item).split('/')
                        if len(a)>10 and len(b)>10:
                            item = round(int(a)/int(b), 4)
                    print(item, end='\t')
                print()
            else:
                item = matrix[i]
                if '/' in str(item):
                    a, b = str(item).split('/')
                    if len(a)>10 and len(b)>10:
                        item = round(int(a)/int(b), 4)
                print(item, end='\t')
        print()

    # random matrix with given parameters
    def random_matrix(cols: int = randint(1,100), rows: int = randint(1,100), start: int = -42, end: int = 42, useFractions: bool = True, square: bool = False, name: str = "M"):
        if square: cols = rows
        matrix = [[0]*cols for _ in range(rows)]
        for i in range(rows):
            for j in range(cols):
                a = randint(start, end)
                b = randint(start, end) if useFractions else 1
                if b == 0: b = 1
                matrix[i][j] = Fraction(a, b)
        return Matrix(matrix, name)


    # summation of two matrices
    def __add__(self, other):
        if not isinstance(other, Matrix):
            raise ValueError("BOTH PARAMETERS SHOULD BE MATRICES FOR SUMMATION")
        if not(self.rows == other.rows or self.cols == other.cols):
            raise ValueError("BOTH MATRICES SHOULD BE ONE DIMENSION FOR SUMMATION")
        if self.rows == 1:
            matrix = [self.matrix[i] + other.matrix[i] for i in range(self.cols)]
        else:
            matrix = [[self.matrix[i][j] + other.matrix[i][j] for j in range(self.cols)] for i in range(self.rows)]
        return Matrix(matrix, f"({self.name}) + ({other.name})")
    
    # subtraction of two matrices
    def __sub__(self, other):
        if not isinstance(other, Matrix):
            raise ValueError("BOTH PARAMETERS SHOULD BE MATRICES FOR SUBTRACTION")
        if not(self.rows == other.rows or self.cols == other.cols):
            raise ValueError("BOTH MATRICES SHOULD BE ONE DIMENSION FOR SUBTRACTION")
        if self.rows == 1:
            matrix = [self.matrix[i] - other.matrix[i] for i in range(self.cols)]
        else:
            matrix = [[self.matrix[i][j] - other.matrix[i][j] for j in range(self.cols)] for i in range(self.rows)]
        return Matrix(matrix, f"({self.name}) - ({other.name})")
    
    # miltiplication of two matrices
    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float) or isinstance(other, Fraction):
            matrix = [[(self.matrix[i][j] if isinstance(self.matrix[i], list) else self.matrix[j])*other for j in range(self.cols)] for i in range(self.rows)]
            if len(matrix) == 1 and isinstance(matrix[0], list): matrix = matrix[0]
            if isinstance(matrix, list):
                return Matrix(matrix, f"({self.name}) * {other}")
            return matrix
        if not isinstance(other, Matrix):
            raise ValueError("BOTH PARAMETERS SHOULD BE MATRICES FOR MULTIPLICATION")
        if not(self.cols == other.rows):
            raise ValueError("COLUMNS FROM FIRST MATRIX SHOULD BE EQUAL TO ROWS FROM SECOND MATRIX FOR MULTIPLICATION")
        matrix = [[0 for _ in range(other.cols)] for _ in range(self.rows)]
        if len(matrix) == 1 and isinstance(matrix[0], list): matrix = matrix[0]
        for i in range(self.rows):
            for j in range(other.cols):
                vec_1 = deepcopy(self.matrix[i] if isinstance(self.matrix[i], list) else self.matrix)
                vec_2 = deepcopy(other.T().matrix[j] if isinstance(other.T().matrix[j], list) else other.T().matrix)
                if isinstance(matrix[i], list):
                    matrix[i][j] = sum([vec_1[k] * vec_2[k] for k in range(self.cols)])
                else:
                    matrix[i] = sum([vec_1[k] * vec_2[k] for k in range(self.cols)])
        if len(matrix) == 1 and isinstance(matrix[0], list): matrix = matrix[0]
        return Matrix(matrix, f"({self.name}) * ({other.name})")
    
    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float) or isinstance(other, Fraction):
            matrix = [[Fraction((self.matrix[i][j] if isinstance(self.matrix[i], list) else self.matrix[j]))/Fraction(other) for j in range(self.cols)] for i in range(self.rows)]
            if len(matrix) == 1 and isinstance(matrix[0], list): matrix = matrix[0]
            if isinstance(matrix, list):
                return Matrix(matrix, f"({self.name}) / {other}")
            else:
                return matrix
        if not isinstance(other, Matrix):
            raise ValueError("BOTH PARAMETERS SHOULD BE MATRICES FOR DIVISION")
        return other.inverse() * self
        

    # text = str or list
    def __init__(self, text, name: str):
        self.name = name
        self.matrix = self.parse_matrix(text) if type(text) is str else list(text)
        if not self.isLegal():
            raise ValueError("MATRIX SHOULD BE NxM")
        self.rows = len(self.matrix) if type(self.matrix[0]) is list else 1
        self.cols = len(self.matrix[0]) if type(self.matrix[0]) is list else len(self.matrix)
