import numpy as np
import math


def splfltlin(pprofile, dx, lc):
    beta = 0.0  # % default by ISO / TS 16610 - 22
    alfa = alfapar(dx, lc)
    n = len(pprofile)
    m1 = np.eye(n)
    q_matrix = qmatlin(n)
    print(q_matrix)

    if beta == 0.0:
        temp = m1 + alfa ** 4 * q_matrix
        temp1 = pprofile.T.conj()
        temp2 = np.linalg.solve(temp, temp1)
        w_profile = temp2.T
    else:
        p_matrix = p_matlin(n)
        temp = m1 + beta * alfa ** 2 * p_matrix
        temp1 = (1 - beta) * alfa ** 4 * q_matrix
        temp2 = temp + temp1
        temp3 = pprofile.T.conj()
        temp4 = np.linalg.solve(temp2, temp3)
        w_profile = temp4.T

    return w_profile


def alfapar(dx, lc):
    if lc <= 0:
        alfa = 0
        raise Exception('ALFAPAR lc is not > 0')
    elif dx <= 0:
        alfa = 0
        raise Exception('ALFAPAR dx is not > 0')
    elif lc < (8 * dx):
        alfa = 0
        raise Exception('ALFAPAR lc is not >= 8*DeltaX')
    else:
        alfa = 1 / (2 * math.sin(dx * math.pi / lc))
    return alfa


def p_matlin(n):
    if n <= 2:
        print('QMATLIN n is not > 4')
        return
    else:
        p_matrix = np.zeros((n, n), dtype=int)

        for i in range(n):
            p_matrix[i, i] = 2
        for i in range(0, n - 1):
            p_matrix[i, i + 1] = -1
            p_matrix[i + 1, i] = -1
        p_matrix[0, 0] = 1
        p_matrix[n - 1, n - 1] = 1
    return p_matrix


def qmatlin(n):
    if n < 4:
        print('QMATLIN n is not > 4')
        return None
    else:
        q_matrix = np.zeros((n, n), dtype=int)
        q_matrix[0, 0] = 1
        q_matrix[0, 1] = -2
        q_matrix[1, 0] = -2
        q_matrix[1, 1] = 5
        q_matrix[n - 1, n - 1] = 1
        q_matrix[n - 1, n - 2] = -2
        q_matrix[n - 2, n - 1] = -2
        q_matrix[n - 2, n - 2] = 5
        for i in range(n - 2):
            q_matrix[i, i + 2] = 1
            q_matrix[i + 2, i] = 1
        for i in range(1, n - 2):
            q_matrix[i, i + 1] = -4
            q_matrix[i + 1, i] = -4
        for i in range(2, n - 2):
            q_matrix[i, i] = 6
    return q_matrix


if __name__ == '__main__':
    m = qmatlin(10)
    p = p_matlin(10)
    # print(p)
    x = np.linspace(-5, 5, 1000)
    errors = np.random.normal(0, 1, x.size)
    y = x ** 2 + errors
    y[500] = 10
    print(y)
    ys = splfltlin(y, 0.01, 2.5)
    print(ys)