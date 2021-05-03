# Abro los ficheros con las coordenadas de los carbonos alfa de la proteína 1 y 2
coord1 = open(input("Ingresar nombre del archivo que contiene las coordenadas de los carbonos alfa de la proteína 1 \n(el archivo debe estar en el mismo directorio que este script): "),"r")
coord2 = open(input("Ingresar nombre del archivo que contiene las coordenadas de los carbonos alfa de la proteína 2 \n(el archivo debe estar en el mismo directorio que este script): "),"r")

# Leo todas las líneas del fichero 1
lineas1 = coord1.readlines()

# Bucle para calcular el centro de masa de la proteína 1
i=0
x=0
y=0
z=0

for linea in lineas1:
    i = i+1

    coordX = float(linea[31:38]) 
    coordY = float(linea[39:46])
    coordZ = float(linea[47:54])

    x = x + coordX
    y = y + coordY
    z = z + coordZ

cmX1 = x/i
cmY1 = y/i
cmZ1 = z/i

# Imprimo las coordenadas del centro de masa de la proteína 1
print("Centro de masa de la proteína 1:",cmX1,cmY1,cmZ1)

# Cierro el fichero de coordenadas de la proteína 1  
coord1.close()

# Leo todas las líneas del fichero 2
lineas2 = coord2.readlines()

# Bucle para calcular el centro de masa de la proteína 2
i=0
x=0
y=0
z=0

for linea in lineas2:
    i = i+1
    
    coordX = float(linea[31:38]) 
    coordY = float(linea[39:46])
    coordZ = float(linea[47:54])

    x = x + coordX
    y = y + coordY
    z = z + coordZ

cmX2 = x/i
cmY2 = y/i
cmZ2 = z/i  

# Imprimo las coordenadas del centro de masa de la proteína 2
print("Centro de masa de la proteína 2:",cmX2,cmY2,cmZ2)

# Cierro el fichero de coordenadas de la proteína 1  
coord2.close()
