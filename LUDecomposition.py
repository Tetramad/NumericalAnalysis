import scipy as sp
import scipy.linalg as linalg


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


def AssertLU(matrix: sp.ndarray, msg=''):
    print(msg)

    L, U = LUDecomposition(matrix)
    if sp.allclose(matrix, L@U, rtol=0):
        print('PASS: LU decomposition')
    else:
        print('FAIL: LU decomposition')

    P, L, U = LUPPDecomposition(matrix)
    if sp.allclose(matrix, P@L@U, rtol=0):
        print('PASS: LUPP decomposition')
    else:
        print('FAIL: LUPP decomposition')

    P, L, U, Q = LUCPDecomposition(matrix)
    if sp.allclose(matrix, P@L@U@Q, rtol=0):
        print('PASS: LUCP decomposition')
    else:
        print('FAIL: LUCP decomposition')

    if (matrix.shape[0] == matrix.shape[1]):
        L = CholeskyDecomposition(matrix)
        if sp.allclose(matrix, L@L.transpose(), rtol=0):
            print('PASS: Cholesky decomposition')
        else:
            print('FAIL: Cholesky decomposition')

    print()


def main():
    A = sp.array([[1, 1, 1],
                  [1, 2, 4],
                  [1, 3, 9]])
    AssertLU(A, 'Low rank Vandermonde matrix')

    A = sp.array([[-1]])
    AssertLU(A, 'Single element matrix')

    A = sp.array([[3, 2, 4],
                  [2, 4, 3]])
    AssertLU(A, 'Non-square matrix')

    A = sp.array([[2, 4],
                  [3, 3],
                  [4, 2]])
    AssertLU(A, 'Non-square matrix')

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
    AssertLU(A, 'High rank Vandermonde matrix')

    A = sp.array([[1, 1, 1],
                  [2, 2, 2],
                  [3, 3, 3]])
    AssertLU(A, 'Singular matrix')


if __name__ == '__main__':
    main()
