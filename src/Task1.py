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
M = 666

h1 = l1 / N1
h2 = l2 / N2
tau = T / M

xi = [i * h1 for i in range(1, N1 + 1)]
yi = [i * h2 for i in range(1, N2 + 1)]
ti = [i * tau for i in range(1, M + 1)]


def getMu(x, y, t):
    return 1


def getF(x, y, t):
    return 0


def getAnalitic():
    global N1, N2, M
    y = np.zeros((M, N1, N2))
    for t in range(M):
        for m in range(N1):
            for n in range(N2):
                y[t][m][n] = exp(xi[m] + yi[n] + 2 * ti[t])
    return y


def getG(x, y, t):
    return exp(xi[x] + yi[y] + 2 * ti[t])


def getTask1Result():
    global N1, N2, M
    y = np.zeros((M, N1, N2))

    for m in range(N1):
        for n in range(N2):
            y[0][m][n] = getG(m, n, 0)
            y[0][m][n] = getG(m, n, 0)

    for t in range(1, M):
        for m in range(N1):
            for n in range(N2):
                y[t][0][n] = getG(0, n, t)
                y[t][N1 - 1][n] = getG(N1 - 1, n, t)
            y[t][m][0] = getG(m, 0, t)
            y[t][m][N2 - 1] = getG(m, N2 - 1, t)

    for j in range(1, M - 1):
        for m in range(1, N1 - 1):
            for n in range(1, N2 - 1):
                y[j][m][n] = tau * (getF(m, n, j) + getMu(m, n, j) *
                                    (((y[j - 1][m - 1][n] - 2 * y[j - 1][m][n] + y[j - 1][m + 1][n]) / (h1 ** 2)) + (
                                        y[j - 1][m][n - 1] - 2 * y[j - 1][m][n] + y[j - 1][m][n + 1]) / (h2 ** 2))) + \
                             y[j - 1][m][n]
    return y


y = getTask1Result()
yA = getAnalitic()
fig1 = pylab.figure()
axes = fig1.gca(projection='3d')
xi, yi = np.meshgrid(xi, yi)
print(xi.shape)
axes.plot_surface(xi, yi, y[ti.index(tau * (M - 1))])
print(np.linalg.norm(y - yA) / np.linalg.norm(yA))
pylab.show()
