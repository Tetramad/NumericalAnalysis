import scipy as sp
import scipy.linalg as linalg
from timeit import timeit


def DoolittleLU(matrix: sp.ndarray, dtype=float):
    n, m = matrix.shape
    l = sp.zeros((n, n), dtype=dtype)
    u = sp.zeros((n, m), dtype=dtype)

    for i in range(n):
        for k in range(i, m):
            s = 0
            for j in range(i):
                s += l[i][j] * u[j][k]

            u[i][k] = matrix[i][k] - s

        l[i][i] = 1
        for k in range(i+1, n):
            s = 0
            for j in range(i):
                s += l[k][j] * u[j][i]
            l[k][i] = (matrix[k][i] - s) / u[i][i]

    return l, u


def CroutLU(matrix: sp.ndarray, dtype=float):
    n, m = matrix.shape
    l = sp.zeros((n, m))
    u = sp.zeros((m, m))

    for j in range(m):
        for i in range(j, n):
            s = 0
            for k in range(j):
                s += l[i][k] * u[k][j]
            l[i][j] = matrix[i][j] - s

        u[j][j] = 1
        for i in range(j+1, m):
            s = 0
            for k in range(j):
                s += l[j][k] * u[k][i]
            u[j][i] = (matrix[j][i] - s) / l[j][j]

    return l, u


def LUDecomposition(matrix: sp.ndarray, dtype=float):
    n, m = matrix.shape
    L = sp.identity(n, dtype)
    U = sp.array(matrix, dtype=dtype)

    for i in range(1, n):
        for j in range(i):
            L[i][j], U[i][j] = U[i][j]/U[j][j], dtype(0)
            for k in range(j+1, m):
                U[i][k] -= U[j][k]*L[i][j]

    return L, U


def LUPPDecomposition(matrix: sp.ndarray, dtype=float):
    n, _ = matrix.shape
    absmatrix = sp.absolute(matrix)
    scales = sp.amax(absmatrix, axis=1)
    p = sp.arange(n)
    p = sorted(p, key=lambda i: absmatrix[i][0]/scales[i], reverse=True)

    P = sp.zeros((n, n))
    for idx, pivot in enumerate(p):
        P[idx][pivot] = 1

    L, U = LUDecomposition(P@matrix, dtype)
    return P.transpose(), L, U


def LUCPDecomposition(matrix: sp.ndarray, dtype=float):
    n, m = matrix.shape
    absmatrix = sp.absolute(matrix)
    rowscales = sp.amax(absmatrix, axis=1)
    colscales = sp.amax(absmatrix, axis=0)
    p = sp.arange(n)
    p = sorted(p, key=lambda i: absmatrix[i][0]/rowscales[i], reverse=True)
    q = sp.arange(m)
    q = sorted(q, key=lambda j: absmatrix[0][j]/colscales[j], reverse=True)

    P = sp.zeros((n, n))
    for idx, pivot in enumerate(p):
        P[idx][pivot] = 1
    Q = sp.zeros((m, m))
    for idx, pivot in enumerate(q):
        Q[idx][pivot] = 1

    L, U = LUDecomposition(P@matrix@Q, dtype)

    return P.transpose(), L, U, Q.transpose()


def CholeskyDecomposition(matrix: sp.ndarray, dtype=float):
    n, _ = matrix.shape
    L = sp.zeros((n, n))

    for i in range(n):
        for j in range(i):
            L[i][j] = matrix[i][j]
            for k in range(j):
                L[i][j] -= L[i][k] * L[j][k]
            L[i][j] /= L[j][j]
        L[i][i] = matrix[i][i]
        for k in range(i):
            L[i][i] -= L[i][k] ** 2
        L[i][i] **= 0.5

    return L


def LDLTDecomposition(matrix: sp.ndarray, dtype=float):
    L, U = LUDecomposition(matrix, dtype)

    return L, sp.diag(sp.diag(U))


def AssertLU(matrix: sp.ndarray, msg=''):
    print(msg)
    P, L, U = linalg.lu(matrix)
    if sp.allclose(matrix, P@L@U, rtol=0):
        print('PASS: linalg.lu')
        print(
            f'time: {timeit(lambda: linalg.lu(matrix), number=1000) * 1000: 3.0f}ms')
    else:
        print('FAIL: linalg.lu')

    L, U = DoolittleLU(matrix)
    if sp.allclose(matrix, L@U, rtol=0):
        print('PASS: DoolittleLU')
        print(
            f'time: {timeit(lambda: DoolittleLU(matrix), number=1000) * 1000: 3.0f}ms')
    else:
        print('FAIL: DoolittleLU')

    L, U = CroutLU(matrix)
    if sp.allclose(matrix, L@U, rtol=0):
        print('PASS: CroutLU')
        print(
            f'time: {timeit(lambda: CroutLU(matrix), number=1000) * 1000: 3.0f}ms')
    else:
        print('FAIL: CroutLU')

    L, U = LUDecomposition(matrix)
    if sp.allclose(matrix, L@U, rtol=0):
        print('PASS: LU decomposition')
        print(
            f'time: {timeit(lambda: LUDecomposition(matrix), number=1000) * 1000: 3.0f}ms')
    else:
        print('FAIL: LU decomposition')

    P, L, U = LUPPDecomposition(matrix)
    if sp.allclose(matrix, P@L@U, rtol=0):
        print('PASS: LUPP decomposition')
        print(
            f'time: {timeit(lambda: LUPPDecomposition(matrix), number=1000) * 1000: 3.0f}ms')
    else:
        print('FAIL: LUPP decomposition')

    P, L, U, Q = LUCPDecomposition(matrix)
    if sp.allclose(matrix, P@L@U@Q, rtol=0):
        print('PASS: LUCP decomposition')
        print(
            f'time: {timeit(lambda: LUCPDecomposition(matrix), number=1000) * 1000: 3.0f}ms')
    else:
        print('FAIL: LUCP decomposition')

    # Symmetric matrix only
    if matrix.shape[0] == matrix.shape[1] and sp.allclose(matrix, sp.transpose(matrix), rtol=0):
        L, D = LDLTDecomposition(matrix)
        if sp.allclose(matrix, L@D@L.transpose(), rtol=0):
            print('PASS: LDLT decomposition')
            print(
                f'time: {timeit(lambda: LDLTDecomposition(matrix), number=1000) * 1000: 3.0f}ms')
        else:
            print('FAIL: LDLT decomposition')

        L = CholeskyDecomposition(matrix)
        if sp.allclose(matrix, L@L.transpose(), rtol=0):
            print('PASS: Cholesky decomposition')
            print(
                f'time: {timeit(lambda: CholeskyDecomposition(matrix), number=1000) * 1000: 3.0f}ms')
        else:
            print('FAIL: Cholesky decomposition')

    print()


def main():
    A = sp.array([[1, 1, 1],
                  [1, 2, 4],
                  [1, 3, 9]])
    AssertLU(A, 'Low rank(3) Vandermonde matrix')

    A = sp.array([[-1]])
    AssertLU(A, 'Single element matrix')

    A = sp.array([[3, 2, 4],
                  [2, 4, 3]])
    AssertLU(A, 'Non-square(2, 3) matrix')

    A = sp.array([[2, 4],
                  [3, 3],
                  [4, 2]])
    AssertLU(A, 'Non-square(3, 2) matrix')

    A = sp.array([[1, -2, 0],
                  [-2, 1, -2],
                  [0, -2, 1]])
    AssertLU(A, 'Symmetric matrix')

    A = sp.array([[4, 12, -16],
                  [12, 37, -43],
                  [-16, -43, 98]])
    AssertLU(A, 'Positive-definite matrix')

    A = sp.array([[1, 2, 3],
                  [2, 4, 7],
                  [3, 3, 3]])
    AssertLU(A, 'Bad condition of naive LU decomposition')

    A = sp.array([[1, 9, 3],
                  [2, 2, 7],
                  [3, 3, 3]])
    AssertLU(A, 'Bad condition of LUPP decomposition')

    A = sp.zeros((16, 16))
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            A[i][j] = (i+1)**j
    AssertLU(A, 'High rank(16) Vandermonde matrix')

    A = sp.zeros((17, 17))
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            A[i][j] = (i+1)**j
    AssertLU(A, 'High rank(17) Vandermode matrix')

    A = sp.array([[1, 1, 1],
                  [2, 2, 2],
                  [3, 3, 3]])
    AssertLU(A, 'Singular matrix')

    A = sp.rand(50, 50)
    AssertLU(A, 'Big(50, 50) matrix')


if __name__ == '__main__':
    main()
