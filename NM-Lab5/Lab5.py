import numpy as np
import sympy as sp
import math
import matplotlib.pyplot as plt

xp = sp.Array([-1.5, -1, -0.5, 0, 0.5, 1, 1.5])


def func_np(x):
    y = 8 * np.sin(1 / np.tan(3)) + (np.sin(5 * x)) ** 2 / (5 * np.cos(10 * x))
    return y


def func(x):
    y = 8 * sp.sin(1 / sp.tan(3)) + (sp.sin(5 * x)) ** 2 / (5 * sp.cos(10 * x))
    return y


def func_six():
    x = sp.Symbol('x')
    r = sp.diff(func(x), x, xp.shape[0])
    return sp.lambdify(x, r, 'numpy')


def newton_up():
    n = xp.shape[0]
    x = sp.Symbol('x')
    res = 0
    div = func(x)
    for i in range(0, n, 1):
        temp = div.evalf(subs={x: xp[i]})
        for j in range(0, i, 1):
            temp *= (x - xp[j])
        res += temp
        div = (div - div.evalf(subs={x: xp[i]})) / (x - xp[i])
    f = open('newton_up.txt', 'w')
    f.write(str(res) + '\n\n\n')
    f.close()
    return sp.lambdify(x, res, 'numpy')


def newton_down():
    n = xp.shape[0]
    x = sp.Symbol('x')
    res = 0
    div = func(x)
    for i in range(n-1, -1, -1):
        temp = div.evalf(subs={x: xp[i]})
        for j in range(n-1, i, -1):
            temp *= (x - xp[j])
        res += temp
        div = (div - div.evalf(subs={x: xp[i]})) / (x - xp[i])
    f = open('newton_down.txt', 'w')
    f.write(str(res) + '\n\n\n')
    f.close()
    return sp.lambdify(x, res, 'numpy')


def lagrange():
    n = xp.shape[0]
    x = sp.Symbol('x')
    pol = 1
    for i in range(0, n, 1):
        pol *= (x - xp[i])
    res = 0
    for i in range(0, n, 1):
        w = pol/(x - xp[i])
        res += func(xp[i])*w/w.evalf(subs={x: xp[i]})
    f = open('lagrange.txt', 'w')
    f.write(str(res) + '\n\n\n')
    f.close()
    return sp.lambdify(x, res, 'numpy')


def spline():
    n = xp.shape[0]
    y = np.zeros(n)
    for i in range(0, n, 1):
        y[i] = func(xp[i])
    a = np.copy(y)
    b = np.zeros(n)
    c = np.zeros(n)
    d = np.zeros(n)
    alpha = np.zeros(n)
    beta = np.zeros(n)
    A = B = C = D = 0

    for i in range(1, n - 1, 1):
        h_i = xp[i] - xp[i-1]
        h_i1 = xp[i+1] - xp[i]

        A = h_i
        C = 2*(h_i + h_i1)
        B = h_i1
        D = 6*((y[i+1] - y[i])/h_i1 - (y[i] - y[i-1])/h_i)
        z = (A * alpha[i-1] + C)
        alpha[i] = -B/z
        beta[i] = (D - A*beta[i-1])/z

    c[n-1] = (D - A*beta[n-2])/(C + A*alpha[n-2])

    for i in range(n-2, 0, -1):
        c[i] = alpha[i] * c[i+1] + beta[i]

    for i in range(n-1, 0, -1):
        h_i = xp[i] - xp[i-1]
        d[i] = (c[i] - c[i-1])/h_i
        b[i] = h_i*(2*c[i] + c[i-1])/6 + (y[i] - y[i-1])/h_i

    f = open('spline.txt', 'w')
    x = sp.Symbol('x')
    res = sp.zeros(n)
    res[0] = a[0] + b[0] * (x - xp[0]) + c[0] / 2 * (x - xp[0]) ** 2 + d[0] / 6 * (x - xp[0]) ** 3
    f.write('If x < ' + str('%.g' % (xp[0])) + ':\n')
    f.write('\t' + str(res[0]) + '\n')

    for i in range(1, n-1, 1):
        res[i] = a[i] + b[i]*(x-xp[i]) + c[i]/2*(x-xp[i])**2 + d[i]/6*(x-xp[i])**3
        f.write('If ' + str('%.g' % (xp[i-1])) + ' < x < ' + str('%.g' % (xp[i])) + ':\n')
        f.write('\t' + str(res[i]) + '\n')

    res[n-1] = a[n-1] + b[n-1] * (x - xp[n-1]) + c[n-1] / 2 * (x - xp[n-1]) ** 2 + d[n-1] / 6 * (x - xp[n-1]) ** 3
    f.write('If x > ' + str('%.g' % (xp[n-2])) + ':\n')
    f.write('\t' + str(res[n-1]) + '\n\n\n')
    f.close()
    return res


def spline_calc(method, arg):
    x = sp.Symbol('x')
    n = method.shape[0]
    i = 0
    while arg > xp[i] and i < n-1:
        i += 1
    f = sp.lambdify(x, method[i], 'numpy')
    z = f(arg)
    return z


def majorant():
    n = xp.shape[0]
    x = sp.Symbol('x')
    w = 1
    for i in range(0, n, 1):
        w *= (x - xp[i])
    w = abs(w)/math.factorial(6)

    ab = np.linspace(np.float(xp[0]), np.float(xp[n-1]))
    der = func_six()
    temp = np.max(abs(der(ab)))

    w *= temp
    f = open('majorant.txt', 'w')
    f.write(str(w) + '\n\n\n')
    f.close()
    return sp.lambdify(x, w, 'numpy')


def draw():
    new_up = newton_up()
    new_down = newton_down()
    lagr = lagrange()
    splin = spline()
    maj = majorant()

    n = xp.shape[0]
    space = np.linspace(np.float(xp[0]), np.float(xp[n-1]), np.int((xp[n-1] - xp[0])*100 + 1))

    f = open('newton_up.txt', 'a')
    f.write('N\tX\t\tf(x)\t\t\t\tNewton(up)\t\t\teps\n')
    for i in range(0, space.shape[0], 10):
        f.write(str(i) + '\t' + str('%.2f' % (space[i])) + '\t' + str(func_np(space[i])) + '\t')
        f.write(str(new_up(space[i])) + '\t' + str(abs(new_up(space[i]) - func_np(space[i]))) + '\n')
    f.close()

    f = open('newton_down.txt', 'a')
    f.write('N\tX\t\tf(x)\t\t\t\tNewton(down)\t\t\teps\n')
    for i in range(0, space.shape[0], 10):
        f.write(str(i) + '\t' + str('%.2f' % (space[i])) + '\t' + str(func_np(space[i])) + '\t')
        f.write(str(new_down(space[i])) + '\t' + str(abs(new_down(space[i]) - func_np(space[i]))) + '\n')
    f.close()

    f = open('lagrange.txt', 'a')
    f.write('N\tX\t\tf(x)\t\t\t\tLagrange\t\t\teps\n')
    for i in range(0, space.shape[0], 10):
        f.write(str(i) + '\t' + str('%.2f' % (space[i])) + '\t' + str(func_np(space[i])) + '\t')
        f.write(str(lagr(space[i])) + '\t' + str(abs(lagr(space[i]) - func_np(space[i]))) + '\n')
    f.close()

    f = open('spline.txt', 'a')
    f.write('N\tX\t\tf(x)\t\t\t\tSpline\t\t\t\teps\n')
    for i in range(0, space.shape[0], 10):
        f.write(str(i) + '\t' + str('%.2f' % (space[i])) + '\t' + str(func_np(space[i])) + '\t')
        f.write(str(spline_calc(splin, space[i])) + '\t' + str(abs(spline_calc(splin, space[i]) - func_np(space[i]))) + '\n')
    f.close()

    f = open('majorant.txt', 'a')
    f.write('N\tX\t\tMajorant\n')
    for i in range(0, space.shape[0], 10):
        f.write(str(i) + '\t' + str('%.2f' % (space[i])) + '\t' + str(maj(space[i])) + '\n')
    f.close()

    plt.plot(space, abs(new_up(space) - func_np(space)), label='Newton (up)')
    plt.plot(space, abs(new_down(space) - func_np(space)), label='Newton (down)')
    plt.plot(space, abs(lagr(space) - func_np(space)), label='Lagrange')
    temp = np.zeros(space.shape[0])
    for i in range(0, space.shape[0], 1):
        temp[i] = spline_calc(splin, space[i])
    plt.plot(space, abs(temp - func_np(space)), label='Splines')
    plt.legend(loc='upper center')
    plt.title('Error')
    plt.grid(True)

    plt.figure()
    plt.plot(space, abs(new_up(space) - func_np(space)), label='Newton (up)')
    plt.plot(space, abs(new_down(space) - func_np(space)), label='Newton (down)')
    plt.plot(space, abs(lagr(space) - func_np(space)), label='Lagrange')
    temp = np.zeros(space.shape[0])
    for i in range(0, space.shape[0], 1):
        temp[i] = spline_calc(splin, space[i])
    plt.plot(space, abs(temp - func_np(space)), label='Splines')
    plt.plot(space, maj(space), label='Majorant')
    plt.legend(loc='upper center')
    plt.title('Error + majorant')
    plt.grid(True)

    plt.show()


draw()
