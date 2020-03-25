#programa que calcula la cantidad de nucleotidos en una secuencia de ADN, el porcentaje de cada uno, el porcentaje GC y el
#total de nucleotidos.

#llamo 'ñ' al archivo que contiene la secuencia nucleotidica.

ñ=open(input('Nombre del archivo .txt: '), 'r')
secuencia=ñ.readlines()

#a las variables 'a', 'c', 'g', 't' y 'u' les asigno el valor de la cantidad de veces que aparece dicho nucleotido en la
#secuencia. llamo 'gc' a la suma de nucleotidos g y c. llamo 'y' al total de nucleotidos.

for nucleotido in secuencia:
    a=nucleotido.count('A')
    c=nucleotido.count('C')
    g=nucleotido.count('G')
    t=nucleotido.count('T')
    u=nucleotido.count('U')
    gc=g+c
    y=a+c+g+t+u
    print('a:',a,' %a:',a/y*100)
    print('c:',c,' %c:',c/y*100)
    print('g:',g,' %g:',g/y*100)
    print('t:',t,' %t:',t/y*100)
    print('u:',u,' %u:',u/y*100)
    print('GC:',gc,'%GC:',gc/y*100)
    print('total:',y)

ñ.close()