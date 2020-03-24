def calcular_media(lista_de_variables):
	total = 0
	for i in lista_de_variables:
		total += i
	resultado = total / len(lista_de_variables)
	return resultado



numero_de_variables = int(input("Ingrese la cantidad de variables con las que desea operar: "))

lista_de_variables = []

for i in range(numero_de_variables):
	lista_de_variables.append(float(input("Ingrese la variable: ")))
	i += i


print()
print(f"Sus variables son: {lista_de_variables}")
print()

media = calcular_media(lista_de_variables)
print(f"La media es: {media}")
