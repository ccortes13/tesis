from Bio import SeqIO
from Bio.Seq import Seq

def adn_puntos(archivo):



    C = (0,1)
    G = (1,1)
    A = (0,0)
    T = (1,0)

    x_0 = (0.5,0.5)


    fa = open(archivo,'r')
    fread = SeqIO.parse(fa, "fasta")

    secuencia = ""

    for record in fread:
        if len(secuencia) < 1500000: 
            secuencia += record.seq
        else:
            break



    def punto_medio(x_0,nucleotido):
        x = (x_0[0] + nucleotido[0])/2.0
        y = (x_0[1] + nucleotido[1])/2.0
        return (x,y)

    def siguiente_punto(x_0,nucleotido):
        if nucleotido == "A":
            return punto_medio(x_0,A)
        elif nucleotido == "C":
            return punto_medio(x_0,C)
        elif nucleotido == "G":
            return punto_medio(x_0,G)
        elif nucleotido == "T":
            return punto_medio(x_0,T)
        else:
            return x_0

    secuencia_puntos = []

    #lineas de prueba
    secuencia = secuencia[:1500000]


    for nuc in secuencia:
        secuencia_puntos.append((x_0[0],x_0[1]))
        x_0 = siguiente_punto(x_0,nuc)

    return secuencia_puntos
  
