def func(x):
    return x**4 - 3*x**3 + x**2 - 2*x - 8


def deriv(x):
    return 4*x**3 - 9*x**2 + 2*x - 2


def half_divide(a, b, f):
    file = open('results.txt','a')
    file.write("Half divide method\n")
    file.write("N\tA\t\tB\t\t|B-A|\t\tX\t\tf(X)\n")
    count = 0
    x = (a+b)/2
    while abs(b-a) >= e and abs(f(x)) >= e:
        x = (a+b)/2
        count += 1
        file.write("%i\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n" % (count, a, b, abs(b-a), x, f(x)))
        a, b = (a, x) if f(a) * f(x) < 0 else (x, b)
    file.write("Number of iterations: " + str(count))
    file.write("\nAnswer: " + str(x) + "\n\n")
    file.close()
    return x


def chord(a, b, f):
    file = open('results.txt', 'a')
    file.write("Chord method\n")
    file.write("N\tA\t\tB\t\t|B-A|\t\tX\t\tf(X)\n")
    count = 0
    x = (a*f(b) - b*f(a))/(f(b) - f(a))
    while abs(b-a) >= e and abs(f(x)) >= e:
        x = (a * f(b) - b * f(a)) / (f(b) - f(a))
        count += 1
        file.write("%i\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n" % (count, a, b, abs(b-a), x, f(x)))
        a, b = (a, x) if f(a) * f(x) < 0 else (x, b)
    file.write("%i\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n" % (count+1, a, b, abs(b-a), x, f(x)))
    file.write("Number of iterations: " + str(count))
    file.write("\nAnswer: " + str(x) + "\n\n")
    file.close()
    return x


def newton(a, f, df):
    file = open('results.txt', 'a')
    file.write("Newton\'s method\n")
    file.write("N\tX(k)\t\tX(k+1)\t\t|X(k)-X(k+1)|\tf(X(k))\n")
    count = 1
    x = a
    y = x - f(x)/df(x)
    file.write("%i\t%.6f\t%.6f\t%.6f\t%.6f\n" % (count, x, y, abs(x-y), f(x)))
    while abs(y-x) >= e and abs(f(x)) >= e:
        count += 1
        x = y
        y = x - f(x) / df(x)
        file.write("%i\t%.6f\t%.6f\t%.6f\t%.6f\n" % (count, x, y, abs(x - y), f(x)))
    file.write("Number of iterations:  " + str(count))
    file.write("\nAnswer: " + str(x) + "\n\n")
    file.close()
    return x


g = open("results.txt", "w")
g.close()

e = 0.00001

l, r = 3, 3.5
x1 = half_divide(l, r, func)
x2 = chord(l, r, func)
x3 = newton(r, func, deriv)
print(x1, x2, x3)

l,r = -1.2, -0.74
x1 = half_divide(l, r, func)
x2 = chord(l, r, func)
x3 = newton(l, func, deriv)
print(x1, x2, x3)
