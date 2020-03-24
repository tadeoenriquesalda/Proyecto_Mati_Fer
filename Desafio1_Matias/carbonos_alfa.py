# Este script abre el pdb que se indique y copia en un txt nuevo solo las lineas de los carbonos alfa.
f = open(input('Nombre del archivo pdb: '), 'r')
g = open(input('Nombre del nuevo archivo: '), 'w')

lineas = f.readlines()

for linea in lineas:
	condition1 = linea[0:4]
	condition2 = linea[13:15]
	if condition1 == 'ATOM':
		if condition2 == 'CA':
			g.write(linea)

g.close()
f.close()
