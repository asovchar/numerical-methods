import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

e = 0.00001
alpha = 0.856
beta = -1.931
a = 0.564
b = -1.362
c = 0.995
d = 0.051


def f1(x, y):
    z = np.cos(x + alpha) + b * y - c
    return z


def f2(x, y):
    z = x + np.sin(y + beta) - d
    return z


def f3(x, y):
    z = np.sin(x + y) + c * x - d
    return z


def f3_der(x, y):
    z = sp.sin(x + y) + c * x - d
    return z


def f4(x, y):
    z = x ** 2 + y ** 2 - 1
   # z = sp.cos(x + y) + c * x - d
    return z


def draw():
    t = np.linspace(-5, 5)
    x, y = np.meshgrid(t, t)

    plt.subplot(1, 2, 1)
    plt.gca().set_aspect('equal')
    plt.grid(True)
    plt.title('System 1')
    plt.contour(x, y, f1(x, y), [0], colors='red')
    plt.contour(x, y, f2(x, y), [0], colors='blue')

    plt.subplot(1, 2, 2)
    plt.gca().set_aspect('equal')
    plt.grid(True)
    plt.title('System 2')
    plt.contour(x, y, f3(x, y), [0], colors='orange')
    plt.contour(x, y, f4(x, y), [0], colors='green')

    plt.show()


def xy_new(x, y):
    p = (c - np.cos(x + alpha)) / b
    q = d - np.sin(y + beta)
    return q, p


def simple_iteration(left, right):
    f = open('simple.txt', 'w')
    count = 0
    old = np.array([left - 1, right + 1])
    new = np.array([left, right])
    answer = np.ones(2)
    while np.linalg.norm(new - old) > e and np.linalg.norm(answer) > e:
        old = np.copy(new)
        new = xy_new(old[0], old[1])
        answer = f1(old[0], old[1]), f2(old[0], old[1])
        count += 1
        f.write('Iteration ' + str(count) + ':\n')
        f.write('\tThe (x, y) vector is ' + str(old) + '\n')
        f.write('\tThe discrepancy vector is ' + str(answer) + '\n')
        f.write('\tThe F(x,y) vector norm is ' + str(np.linalg.norm(answer)) + '\n\n')
    f.write('The answer is ' + str(new))
    f.close()
    return new


def df3x(p, q):
    x = sp.Symbol('x')
    y = sp.Symbol('y')
    r = sp.diff(f3_der(x, y), x)
    return r.evalf(subs={x: p, y: q})


def df3y(p, q):
    x = sp.Symbol('x')
    y = sp.Symbol('y')
    r = sp.diff(f3_der(x, y), y)
    return r.evalf(subs={x: p, y: q})


def df4x(p, q):
    x = sp.Symbol('x')
    y = sp.Symbol('y')
    r = sp.diff(f4(x, y), x)
    return r.evalf(subs={x: p, y: q})


def df4y(p, q):
    x = sp.Symbol('x')
    y = sp.Symbol('y')
    r = sp.diff(f4(x, y), y)
    return r.evalf(subs={x: p, y: q})


def reverse(matr):
    det = matr[0, 0] * matr[1, 1] - matr[0, 1] * matr[1, 0]
    temp = matr[0, 0]
    matr[0, 0] = matr[1, 1]
    matr[1, 1] = temp
    matr[0, 1] *= -1
    matr[1, 0] *= -1
    matr *= 1/det
    return matr


def newton(left, right):
    f = open('newton.txt', 'a')
    count = 0
    old = np.array([left - 1, right + 1])
    new = np.array([left, right])
    answer = np.ones(2)
    deriv = np.array([[np.float(df3x(left, right)), np.float(df3y(left, right))],
                      [np.float(df4x(left, right)), np.float(df4y(left, right))]])
    deriv_rev = reverse(deriv)
    while np.linalg.norm(new - old) > e and np.linalg.norm(answer) > e:
        old = np.copy(new)
        answer = f3(old[0], old[1]), f4(old[0], old[1])
        new = old - np.dot(deriv_rev, answer)
        count += 1
        f.write('Iteration ' + str(count) + ':\n')
        f.write('\tThe (x, y) vector is ' + str(old) + '\n')
        f.write('\tThe discrepancy vector is ' + str(answer) + '\n')
        f.write('\tThe F(x,y) vector norm is ' + str(np.linalg.norm(answer)) + '\n\n')
    f.write('The answer is ' + str(new) + '\n\n************************************************************\n\n')
    f.close()
    return new


g = open('newton.txt', 'w')
g.close()

res1 = simple_iteration(0.6, -0.6)
res2 = newton(-0.4, 0.9)
res3 = newton(0.5, -0.9)
draw()
