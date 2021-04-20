# Este script abre el pdb que se indique y copia en un txt nuevo solo las lineas de los carbonos alfa.

# Abro el fichero con el pdb para su lectura.
f = open(input('Nombre del archivo pdb: '), 'r')
# Abro un nuevo fichero para editarlo con las líneas de los carbonos alfa.
g = open(input('Nombre del nuevo archivo: '), 'w')

# Leo todas las líneas del fichero con el pdb.
lineas = f.readlines()

# Bucle para buscar las líneas de los carbonos alfa y escribirlas en el fichero nuevo.
for linea in lineas:
	condition1 = linea[0:4]
	condition2 = linea[13:15]
	if condition1 == 'ATOM':
		if condition2 == 'CA':
			g.write(linea)

# Cierro los ficheros abiertos.
g.close()
f.close()
