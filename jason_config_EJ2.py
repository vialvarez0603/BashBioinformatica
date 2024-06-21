#!/usr/bin/env python3

import json
import sys

# Verifica si se proporcionó un argumento para path
if len(sys.argv) < 2:
    print("Error: Debe proporcionar el nombre del archivo JSON como argumento.")
    sys.exit(1)

# Obtiene el argumento para path
path = sys.argv[1]

print('Parámetros Default:')
print('\tTipo: "blastp"')
print('\tNúmero de resultados: 500')
print('\tTipo de matriz: BLOSUM80')
print('\tUmbral para el E-value: 10')
print('\tNombre del archivo de salida: "RESULTADO_BLAST.txt"')

choice = int(input('Para usar la configuración default presionar "0"; Si desea ajustar los parámetros, presione "1":\n'))
while choice != 0 and choice != 1:
    choice = int(input('Vuelva a intentar. Ingrese "0" o "1": '))

if choice == 0:
    tipo = "blastp"
    resultados = 500
    matriz = "BLOSUM80"
    e_thres = 10
    output_file = "RESULTADO_BLAST.txt"
else:
    # Define los parámetros de configuración
    tipo = input('Tipo: ')
    resultados = int(input('Número de resultados: '))
    matriz = input('Tipo de matriz: ')
    e_thres = float(input('Umbral para el E-value: '))
    output_file = input('Nombre del archivo de salida: ')

config = {
    "tipo": tipo,
    "resultados": resultados,
    "matriz": matriz,
    "e_threshold": e_thres,
    "output_path": output_file
}

# Escribe la configuración en un archivo JSON
with open(path, 'w') as config_file:
    json.dump(config, config_file, indent=4)

print(f"Archivo de configuración '{path}' creado con éxito.")
