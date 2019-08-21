from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt




def dibuja_fractal(archivo,n=1):

    C = (0,n)
    G = (n,n)
    A = (0,0)
    T = (n,0)

    x_0 = (n/2,n/2)

    def siguiente_cuadro(x_0,nuc):
        return ((x_0[0]+nuc[0])/2,(x_0[1]+nuc[1])/2)

    def sig_letra(x_0,nt):
        if nt == "C":
            return siguiente_cuadro(x_0,C)
        elif nt == "G":
            return siguiente_cuadro(x_0,G)
        elif nt == "A":
            return siguiente_cuadro(x_0,A)
        elif nt == "T":
            return siguiente_cuadro(x_0,T)
        else:
            return x_0


    fa = open(archivo,'r')
    fread = SeqIO.parse(fa,"fasta")

    secuencia = ""

    for record in fread:
        secuencia += record.seq

    #jumps = 0.022
    jumps = 1.0/3077

    #ceros o unos dependiendo si es para el dibujo o para la metrica
    fractal = np.ones(n**2)
    fractal = fractal.reshape(n,n)

    for nt in secuencia:
        x_0 = sig_letra(x_0,nt) 
        fractal[x_0[0]][x_0[1]] -= jumps

    #fractal = np.reshape(fractal,n**2)
    #print max(fractal)
    print fractal


    plt.gray()
    plt.imshow(fractal)
    plt.show()
    
dibuja_fractal("simulans.txt",256)
