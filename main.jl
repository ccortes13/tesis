using FastaIO
using Glob
using JuMP
using Clp
using Images, Colors

#Ojo cambio de vertices
n = 256
A = (n,0)
C = (0,n)
G = (n,n)
T = (0,0)

x_0 = (n/2,n/2)

sigCuadro(x_0,nuc) = (ceil((x_0[1]+nuc[1])/2),ceil((x_0[2]+nuc[2])/2))

function sigLetra(x_0,nt)
    if nt == 'A'
        x = A
    elseif nt == 'C'
        x = C
    elseif nt == 'G'
        x = G
    elseif nt == 'T'
        x = T
    else
        x = x_0
    end
    
    return sigCuadro(x_0,x)
end


#Dado un archivo de texto en formato fasta
#regresa toda la secuencia de nucleótidos en un string
function creaSeq(archivo)
    secuencia = ""
    for (name, seq) in FastaReader(archivo)
        secuencia = secuencia*seq
    end
    
    return secuencia
end 

#Regresa Matriz de nxn que representa 
#la imagen fractal en escala de grises

function creaMatriz(secuencia)
    jump = 1/length(secuencia)
    fractal_mat = zeros(n,n)
    x_0 = (n/2,n/2)
    
    for i in secuencia
        x_0 = sigLetra(x_0,i)
        x = Int(x_0[1])
        y = Int(x_0[2])
        fractal_mat[x,y] += jump
    end
    
    return fractal_mat
end

#Toma un archivo fasta y regresa matriz que
#representa al fractal
function archivoMatriz(archivo)
    sec = creaSeq(archivo)
    mat = creaMatriz(sec)
    return mat
end

nombres = glob("*.fsa_nt")
matrices = map(archivoMatriz,nombres)
println(length(matrices))

m = Model(solver = ClpSolver())

@variable(m, x[1:n,1:n] >= 0)

@constraint(m,const1[i=1:n-1,j=1:n], x[i,j]-x[i+1,j] <= 1)
@constraint(m,const2[i=2:n,j=1:n], x[i,j]-x[i-1,j] <= 1)
@constraint(m,const3[i=1:n,j=1:n-1], x[i,j]-x[i,j+1] <= 1)
@constraint(m,const4[i=1:n,j=2:n], x[i,j]-x[i,j-1] <= 1)

function optimiza(C)
    if sum(C) < 0
        C *= -1
    end
    @objective(m, Max, sum(C[i,j]*x[i,j] for i=1:n for j=1:n))
    status = solve(m)
    res = getobjectivevalue(m)
    return res              
end

function renglon(i,matrices)
    test = map(x -> matrices[i]-x,matrices)
    resultados = map(optimiza,test[i+1:length(matrices)])
    return vcat(zeros(i),resultados)
end

for i in 1:8
    println(renglon(i,matrices))
end

for i in 1:8
    println(renglon(i,matrices))
end

print(nombres)
