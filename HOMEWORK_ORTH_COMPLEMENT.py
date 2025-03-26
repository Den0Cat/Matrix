from Matrix import Matrix, runtime
from fractions import Fraction

@runtime
def main():
    # A = Matrix.random_matrix(cols=10, rows=10)     # for tests

    A = Matrix(
        """
        1/2   2   3/4   4   5
        6   7   8   9/5   10
        0   0   0   0   0
        11  12  13  14  15""",
        name="A"
    )

    print(f"{'GIVEN MATRIX':=^50}")
    A.print_matrix()
    print(f"{'':=^50}")

    print(f"{'NULLSPACE OF THE GIVEN MATRIX (N(A) ⟂ C(A^T))':=^50}")
    A.nullspace().print_matrix()
    print(f"{'':=^50}")

    print(f"{'LEFT NULLSPACE OF THE GIVEN MATRIX (N(A^T) ⟂ C(A))':=^50}")
    A.left_nullspace().print_matrix()
    print(f"{'':=^50}")

    print(f"{'ROW SPACE (BASIS) OF THE GIVEN MATRIX (C(A^T) ⟂ N(A))':=^50}")
    RREF, pivots = A.rref()
    pivot_rows = [A.matrix[i] for i in range(len(RREF.matrix)) if any(RREF.matrix[i])]
    Matrix(pivot_rows, "C(A^T)").T().print_matrix()
    print(f"{'':=^50}")
    
    

if __name__ == "__main__":
    try:
        main()
    except TypeError:
        pass
