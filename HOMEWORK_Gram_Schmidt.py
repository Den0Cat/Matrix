from Matrix import Matrix, runtime
from fractions import Fraction

@runtime
def main():
    #A = Matrix.random_matrix(cols=10, rows=10)     # for tests

    A = Matrix(
        """
        0   3   5   7   1
        6   1   15  9   10
        10  3   2   1   15
        2   3   4   5   6
        7   8   9   10  13
        18  33  4   34  29
        """,
        name="A"
    )

    A.print_matrix()

    Q, R = A.Gram_Schimdt()
    Q.print_matrix()
    R.print_matrix()

    (Q*R).print_matrix()
    

if __name__ == "__main__":
    try:
        main()
    except TypeError:
        pass
