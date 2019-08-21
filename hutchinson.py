import numpy as np 


def tableu(n):
    pares = []
    for i in range(n):
        pares += [(i,j) for j in range(n)]

    A = []

    for i in range(n):
        for j in range(n):
            x = pares.index((i,j))
            if i+1 < n:
                a = np.zeros(n*n)
                y = pares.index((i+1,j))
                a[x] = 1
                a[y] = -1 
                A.append(a) 
                A.append(a*-1)
            if j+1 < n:
                a = np.zeros(n*n)
                y = pares.index((i,j+1))
                a[x] = 1
                a[y] = -1
                A.append(a)
                A.append(a*-1)
    return A










