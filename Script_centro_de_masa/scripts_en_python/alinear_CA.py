#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Limpio los pdb y dejo sólo los carbonos alfa
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Abro los ficheros pdbs para su lectura.
proteína1 = open(input('Nombre del archivo pdb proteína 1: '), 'r')
proteína2 = open(input('Nombre del archivo pdb proteína 2: '), 'r')

# Abro un nuevo fichero para cada pdb para editarlo con las líneas de los carbonos alfa.
C_alfa_1 = open("CA_proteína1", 'w')
C_alfa_2 = open("CA_proteína2", 'w')

# Leo todas las líneas del fichero con el pdb.
lineas_1 = proteína1.readlines()


# Bucles para buscar las líneas de los carbonos alfa y escribirlas en el fichero nuevo de cada proteína.
for linea in lineas_1:
	condition1 = linea[0:4]
	condition2 = linea[13:15]
	if condition1 == 'ATOM':
		if condition2 == 'CA':
			C_alfa_1.write(linea)

lineas_2 = proteína2.readlines()

for linea in lineas_2:
	condition1 = linea[0:4]
	condition2 = linea[13:15]
	if condition1 == 'ATOM':
		if condition2 == 'CA':
			C_alfa_2.write(linea)

# Cierro los archivos pdb.
proteína1.close()
proteína2.close()

# Cierro los ficheros de los CA.
C_alfa_1.close()
C_alfa_2.close()

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculo el centro de masa de la proteína 1
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Abro para lectura los ficheros de los CA.
coord1 = open("CA_proteína1", "r")
coord2 = open("CA_proteína2", "r")
coord2_alineadas = open("proteína2_alineada_proteína1", "w")

# Leo todas las líneas del fichero de CA 1
lineas1 = coord1.readlines()

# Bucle para calcular el centro de masa de la proteína 1
i=0
x=0
y=0
z=0

for linea in lineas1:

    coordX = float(linea[31:38]) 
    coordY = float(linea[39:46])
    coordZ = float(linea[47:54])

    x = x + coordX
    y = y + coordY
    z = z + coordZ

    i = i+1

cmX1 = x/i
cmY1 = y/i
cmZ1 = z/i

# Imprimo las coordenadas del centro de masa de la proteína 1
print("Centro de masa de la proteína 1:",cmX1,cmY1,cmZ1)

# Cierro el fichero de coordenadas de la proteína 1  
coord1.close()

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Traslado las coordenadas de los CA de la proteína 2 al centro de masa de la proteína 1
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Leo todas las líneas del fichero de CA 2
lineas2 = coord2.readlines()

# Bucle para trasladar todos los átomos de la proteína 2 al centro de masa de la proteína 1

x=0
y=0
z=0

for linea in lineas2:

    coord2_alineadas.write(linea[0:30])

    coordX = float(linea[31:38]) 
    coordY = float(linea[39:46])
    coordZ = float(linea[47:54])

    x = coordX - cmX1
    y = coordY - cmY1
    z = coordZ - cmZ1

    coord2_alineadas.write(f"{x:.3f}" + " ")
    coord2_alineadas.write(f"{y:.3f}" + " ")
    coord2_alineadas.write(f"{z:.3f}" + " ")
    coord2_alineadas.write("\n")
    


# Cierro el fichero de coordenadas de la proteína 2 y el de la proteína 2 alineada 
coord2.close()
coord2_alineadas.close()
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
