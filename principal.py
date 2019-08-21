import multiprocessing as mp
from scipy.spatial.distance import directed_hausdorff
import caos_secuencia as c
import glob
import numpy as np

def haus(A,B):
    
    if A!=B:
        return (max(directed_hausdorff(A,B)[0],directed_hausdorff(B,A)[0]))*100

    else:
        return 0.0



def main():
    moscas = glob.glob("*.txt")
    print moscas

    sets = [c.adn_puntos(mosca) for mosca in moscas]
    
    #Metodo "naive"
    for i in range(len(sets)-1):
        dist = [haus(sets[i],s) for s in sets[i+1:]]
        print dist


    """
    pool = mp.Pool(processes=4)
    results = pool.map(c.adn_puntos,moscas)
    pool = mp.Pool(processes=4)
    dist = [pool.apply(haus, args=(results[0],r)) for r in results[1:]]
    print dist
    """  
    
     



main()

