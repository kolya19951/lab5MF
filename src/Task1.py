import pylab
import numpy as np
from math import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

l1 = 2
l2 = 2
T = 2

N1 = 20
N2 = 20
M = 1000

h1 = l1 / N1
h2 = l2 / N2
tau = round(T / M, 3)

xi = [i * h1 for i in range(1, N1 + 1)]
yi = [i * h2 for i in range(1, N2 + 1)]
ti = [i * tau for i in range(1, M + 1)]


def getMu(x, y, t):
    return 1


def getF(x, y, t):
    return 0


def tma(a, b, c, d):
    n = len(b)
    alp = np.zeros(n)
    bet = np.zeros(n)
    x = np.zeros(n)
    alp[0] = -c[0] / b[0]
    bet[0] = d[0] / b[0]
    np.insert(a, 0, 0)
    np.append(c, 0)

    for i in range(1, n):
        alp[i] = -(c[i] / (a[i] * alp[i - 1] + b[i]))
        bet[i] = (d[i] - a[i] * bet[i - 1]) / (a[i] * alp[i - 1] + b[i])

    x[n - 1] = (d[n - 1] - a[n - 1] * bet[n - 2]) / (a[n - 1] * alp[n - 2] + b[n - 1])

    for i in range(n - 2, -1, -1):
        x[i] = alp[i] * x[i + 1] + bet[i]

    return x


# #1
def getAnalitic():
    global N1, N2, M
    y = np.zeros((M, N2, N1))
    for t in range(M):
        for m in range(N2):
            for n in range(N1):
                y[t][m][n] = exp(xi[n] + yi[m] + 2 * ti[t])
    return y


# 2
# def getAnalitic():
#     global N1, N2, M
#     y = np.zeros((M, N2, N1))
#     for t in range(M):
#         for m in range(N2):
#             for n in range(N1):
#                 y[t][m][n] = cos(xi[n]) * sinh(yi[m]) * exp(-4 * ti[t])
#     return y

# 3
# def getAnalitic():
#     global N1, N2, M
#     y = np.zeros((M, N2, N1))
#     for t in range(M):
#         for m in range(N2):
#             for n in range(N1):
#                 y[t][m][n] = cos(xi[n] + pi / 2) * cosh(2 * yi[m]) * exp(3 * ti[t])
#     return y


# 1
# def getG(x, y, t):
#     if type(t) is int:
#         return exp(xi[x] + yi[y] + 2 * ti[t])
#     if type(t) is float:
#         return exp(xi[x] + yi[y] + 2 * t)
# #
# #3
def getG(x, y, t):
    if type(t) is int:
        return cos(xi[x] + pi / 2) * cosh(2 * yi[y]) * exp(3 * ti[t])
    if type(t) is float:
        return cos(xi[x] + pi / 2) * cosh(2 * yi[y]) * exp(3 * t)


# 2
# def getG(x, y, t):
#     if type(t) is int:
#         return cos(xi[x]) * sinh(yi[y]) * exp(-4 * ti[t])
#     if type(t) is float:
#         return cos(xi[x]) * sinh(yi[y]) * exp(-4 * t)


def getTask1Result():
    global N1, N2, M
    y = np.zeros((M, N2, N1))

    for m in range(N1):
        for n in range(N2):
            y[0][n][m] = getG(m, n, 0)
            y[0][n][m] = getG(m, n, 0)

    for t in range(1, M):
        for m in range(N1):
            for n in range(N2):
                y[t][n][0] = getG(0, n, t)
                y[t][n][N1 - 1] = getG(N1 - 1, n, t)
            y[t][0][m] = getG(m, 0, t)
            y[t][N2 - 1][m] = getG(m, N2 - 1, t)

    for j in range(1, M):
        for m in range(1, N1 - 1):
            for n in range(1, N2 - 1):
                # y[j][n][m] = tau * (getF(m, n, j) + getMu(m, n, j) *
                #                     (((y[j - 1][n][m - 1] - 2 * y[j - 1][n][m] + y[j - 1][n][m + 1]) / (h1 ** 2)) + (
                #                         y[j - 1][n - 1][m] - 2 * y[j - 1][n][m] + y[j - 1][n + 1][m]) / (h2 ** 2))) + \
                #              y[j - 1][n][m]
                y[j][m][n] = tau * (getF(n, m, j) + getMu(n, m, j) *
                                    (((y[j - 1][m - 1][n] - 2 * y[j - 1][m][n] + y[j - 1][m + 1][n]) / (h1 ** 2)) + (
                                        y[j - 1][m][n - 1] - 2 * y[j - 1][m][n] + y[j - 1][m][n + 1]) / (h2 ** 2))) + \
                             y[j - 1][m][n]
    return y


def LTScheme():
    global N1, N2, M
    y = np.zeros((M, N2, N1))
    print(y.shape)
    for m in range(N1):
        for n in range(N2):
            y[0][n][m] = getG(m, n, 0)
            y[0][n][m] = getG(m, n, 0)

    for t in range(1, M):
        for m in range(N1):
            for n in range(N2):
                y[t][n][0] = getG(0, n, t)
                y[t][n][N1 - 1] = getG(N1 - 1, n, t)
            y[t][0][m] = getG(m, 0, t)
            y[t][N2 - 1][m] = getG(m, N2 - 1, t)

    for j in range(M - 1):
        tj = ti[j] + tau / 2
        YFict = np.zeros((N2, N1))
        for m in range(N1):
            for n in range(N2):
                YFict[n][0] = getG(0, n, tj)
                YFict[n][N1 - 1] = getG(N1 - 1, n, tj)
            YFict[0][m] = getG(m, 0, tj)
            YFict[N2 - 1][m] = getG(m, N2 - 1, tj)
        for n in range(1, N2 - 1):
            a = np.zeros(N2 - 2)
            b = np.zeros(N2 - 2)
            c = np.zeros(N2 - 2)
            d = np.zeros(N2 - 2)
            for m in range(1, N1 - 1):
                KN1 = (tau / (2 * ((h1) ** 2))) * getMu(xi[m], yi[n], tj)
                KN2 = (tau / (2 * ((h2) ** 2))) * getMu(xi[m], yi[n], tj)
                if m == 1:
                    b[m - 1] = 1 + (tau / ((h1) ** 2)) * getMu(m, n, tj)
                    c[m - 1] = -KN1
                    d[m - 1] = KN1 * getG(m - 1, n, tj) + (y[j][n - 1][m] + y[j][n + 1][m]) * KN2 + y[j][n][m] * (
                        1 - 2 * KN2)
                elif m == (N1 - 2):
                    a[m - 1] = -KN1
                    c[len(c) - 2] = -KN2
                    b[m - 1] = 1 + (tau / ((h1) ** 2)) * getMu(m, n, tj)
                    d[m - 1] = KN1 * getG(m + 1, n, tj) + (y[j][n - 1][m] + y[j][n + 1][m]) * KN2 + y[j][n][m] * (
                        1 - 2 * KN2)
                    YFict[n, 1: N2 - 1] = tma(a, b, c, d)
                else:
                    a[m - 1] = -KN1
                    b[m - 1] = 1 + (tau / (h1 ** 2)) * getMu(xi[m], yi[n], tj)
                    c[m - 1] = -KN1
                    d[m - 1] = (y[j][n - 1][m] + y[j][n + 1][m]) * KN2 + y[j][n][m] * (1 - 2 * KN2)

        for m in range(1, N2 - 1):
            a = np.zeros(N1 - 2)
            b = np.zeros(N1 - 2)
            c = np.zeros(N1 - 2)
            d = np.zeros(N1 - 2)
            for n in range(1, N1 - 1):
                KN1 = (tau / (2 * ((h1) ** 2))) * getMu(m, n, j + 1)
                KN2 = (tau / (2 * ((h2) ** 2))) * getMu(m, n, j + 1)
                if n == 1:
                    b[n - 1] = 1 + (tau / (h2 ** 2)) * getMu(m, n, j + 1)
                    c[n - 1] = -KN2
                    d[n - 1] = KN2 * y[j + 1][n - 1][m] + (YFict[n][m - 1] + YFict[n][m + 1]) * KN1 + YFict[n][m] * (
                        1 - 2 * KN1)
                elif n == (N1 - 2):
                    a[n - 1] = -KN2
                    b[n - 1] = 1 + (tau / (h2 ** 2)) * getMu(m, n, j + 1)
                    d[n - 1] = KN2 * y[j + 1][n + 1][m] + (YFict[n][m - 1] + YFict[n][m + 1]) * KN1 + YFict[n][m] * (
                        1 - 2 * KN1)
                    y[j + 1][1:N1 - 1, m] = tma(a, b, c, d)
                else:
                    a[n - 1] = -KN2
                    b[n - 1] = 1 + (tau / (h2 ** 2)) * getMu(m, n, j + 1)
                    c[n - 1] = -KN2
                    d[n - 1] = (YFict[n][m - 1] + YFict[n][m + 1]) * KN1 + YFict[n][m] * (1 - 2 * KN1)
    return y


y = getTask1Result()
# yA = getAnalitic()
fig1 = pylab.figure()
axes = fig1.gca(projection='3d')
xi, yi = np.meshgrid(xi, yi)
print(xi.shape)
axes.plot_surface(xi, yi, y[ti.index(0.998)], cmap='viridis')
# axes.plot_surface(xi, yi, yA[ti.index(1.998)], cmap='viridis')
# print(np.linalg.norm(y - yA) / np.linalg.norm(yA))
pylab.show()
