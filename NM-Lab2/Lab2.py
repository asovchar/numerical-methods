import numpy


def triangle(matr, vect):
    file = open('iterations.txt', 'a')
    file.write('***Creating triangular matrixes***\n\n')
    n, m = matr.shape
    right = numpy.copy(matr)
    left = numpy.eye(n)
    det = 1
    for k in range(0, n, 1):
        #c = 0
        if right[k][k] == 0:
            c = k+1
            if c == n:
                return 'The system can not be solved', ' ', 0
            while right[c][k] == 0:
                c += 1
                if c == n:
                    return 'The system can not be solved', ' ', 0
            for j in range(0, n, 1):
                temp = right[k][j]
                right[k][j] = right[c][j]
                right[c][j] = temp
            temp = vect[c]
            vect[c] = vect[k]
            vect[k] = temp
            det *= -1
        for i in range(k + 1, n, 1):
            left[i][k] = right[i][k] / right[k][k]
            right[i] = right[i] - right[k] / right[k][k] * right[i][k]
        '''if c != 0:
            for i in range(k, n, 1):
                temp = left[k][i]
                left[k][i] = left[c][i]
                left[c][i] = temp'''
        file.write('The ' + str(k+1) + ' iteration\n')
        file.write('Upper matrix:\n' + str(right) + '\n')
        file.write('Lower matrix:\n' + str(left) + '\n\n')
    file.close()
    for i in range(0, n, 1):
        det *= right[i][i]
    return right, left, det


def gauss(matr, vect):
    file = open('iterations.txt', 'a')
    right, left, det = triangle(matr, vect)
    if type(right) is str:
        return right, det
    m, n = matr.shape
    ans = numpy.zeros(n)
    res = numpy.zeros(n)
    file.write("\n***The reverse move of Gauss method***\n\n")
    file.write('The solution of the L*ksi=b equality\n')
    for i in range(0, n, 1):
        ans[i] = (vect[i] - numpy.sum(left[i] * ans)) #/ left[i][i]
        file.write('The ksi vector on the ' + str(i+1) + ' iteration: ' + str(ans) + '\n')
    file.write('\nThe solution of the R*X=ksi equality\n')
    for i in range(n-1, -1, -1):
        res[i] = (ans[i] - numpy.sum(right[i] * res)) / right[i][i]
        file.write('The X vector on the ' + str(n-i) + ' iteration: ' + str(res) + '\n')
    file.write('\nThe final answer is: ' + str(res) + '\n\n\n*********************************************************************************************************************\n\n\n\n\n')
    file.close()
    return res, det


def reverse(matr):
    file = open('iterations.txt', 'a')
    n, m = matr.shape
    one = numpy.eye(n)
    ans = numpy.zeros((n, n))
    for i in range(0, n, 1):
        res, temp = gauss(matr, one[i])
        for j in range(0, n, 1):
            ans[j][i] = res[j]
    file.close()
    return ans


g = open('iterations.txt', 'w')
g.close()

'''A = numpy.array([[6.59, 1.08, 0.99, 1.195, -0.21],
                [1.12, 3.83, 1.3, -1.63, -0.48],
                [0.95, -2.46, 5.77, 2.1, -0.017],
                [1.285, 0.16, 2.1, 5.77, 2],
                [0.69, -0.18, 0.283, -1, 4]])
b = numpy.array([2.1, 0.36, -0.43, 1.84, -0.27])'''
A = numpy.array([[1, 1, 5],
              [2, 2, 6],
              [1, -5, -1]])
b = numpy.array([7, 10, -5])


answer, determinant = gauss(A, b)
rev = reverse(A)
disc = b - A@answer
eye = A@rev

g = open('results.txt', 'w')
if answer is str:
    g.write(str(answer) + '\n\n')
    g.write('The determinant is' + str(determinant) + '\n\n')
else:
    g.write('The answer is ' + str(answer) + '\n\n')
    g.write('The determinant is ' + str(determinant) + '\n\n')
    g.write('The reverse matrix is\n' + str(rev) + '\n\n')
    g.write('The vector of discrepancy is [ ')

    size = disc.shape
    for i in range(0, size[0], 1):
        g.write(str('%.6e' % disc[i]) + '\t')
    g.write(']\n\n')
    g.write('The A * A^-1 matrix is\n[')
    for i in range(0, size[0], 1):
        g.write('[ ')
        for j in range(0, size[0], 1):
            g.write(str('%.6e' % (eye[i][j])) + '\t')
        g.write(']\n')
