#Ojo checar bien el tipo de la letra (ie si en string o char)
#Tambien revisar que vertice corresponde a cada letra, por ahora seguimos deschav99
using FastaIO
using Distances
using Glob


A = (0,0)
C = (0,1)
G = (1,1)
T = (1,0)

puntoMedio(x,y) = [(x[1]+y[1])/2,(x[2]+y[2])/2]

#Las matrices x,y son de tamaño nx2 pero el cálculo se hace con las transpuestas
function haus(x,y)
    D = pairwise(Euclidean(),x',y')
    daB = maximum(minimum(D,2))
    dbA = maximum(minimum(D,1))
    result = max(daB,dbA)
    return result
end



function siguientePunto(x_0,nucleotido)
    if nucleotido == 'A'
        x_1 = A
    elseif nucleotido == 'C'
        x_1 = C
    elseif nucleotido == 'G'
        x_1 = G
    elseif nucleotido == 'T'
        x_1 = T
    else 
        x_1 = x_0
    end 
    
    return puntoMedio(x_0,x_1)
end


#dado un archivo con formato Fasta regresa solo los primeros 1500000 nucleotidos
#OJO aqui con las longitudes de las cadenas
function creaSeq(archivo)
    secuencia = ""
    for (name, seq) in FastaReader(archivo)
        #while length(secuencia) < 100000
        secuencia = secuencia*seq
        #end 
    end
    return secuencia[1:100000]
end 

#Dada una cadena de ADN regresa los puntos en forma de matriz nx2
function adnPuntos(adn)
    x_0 = [0.5,0.5]
    puntos = [x_0[1] x_0[2]]
    for i in adn
        x_0 = siguientePunto(x_0,i)
        puntos = [puntos; [x_0[1] x_0[2]]]
    end
    #regresa una matriz de tamaño nx2
    return puntos
end


function main()
    nombres = glob("*.txt")
    secuencias = []
    for i in nombres
       segmento = creaSeq(i)
       mat = adnPuntos(segmento)
       push!(secuencias,mat)
    end
    l = length(secuencias)
    distancias = zeros(l,l)
    for i in 1:l
        for j in 1:l
	    if (i != j) & (distancias[i,j] == 0)
	        dist = haus(secuencias[i],secuencias[j])
		distancias[i,j] = dist
		distancias[j,i] = dist
	    end
	end
    end

    print(distancias)
end


main()

