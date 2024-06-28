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
print('\tApproach ([1]: Puntuamos gaps y missmatch con -1, y match con +2. [2]: Puntuamos sólo gaps y missmatch con -1.): 2')
print('\tNúmero de secuencias en el alineamiento: 10')
print('\tNombre del archivo de salida: "RESULTADO_ALINEAMIENTO.fasta"')

choice = int(input('Para usar la configuración default presionar "0"; Si desea ajustar los parámetros, presione "1":\n'))
while choice != 0 and choice != 1:
    choice = int(input('Vuelva a intentar. Ingrese "0" o "1": '))

if choice == 0:
    approach = 2
    num = 10
    output_file = "Resultados/Resultados_Ejercicio3/RESULTADO_ALINEAMIENTO.fasta"
else:
    # Define los parámetros de configuración
    approach = int(input('Approach: '))
    while approach not in [1, 2]:
        approach = int(input('Elija un valor apropiado para el approach (1 o 2). Approach: '))
    num = int(input('Número de secuencias en el alineamiento: '))
    output_file = 'Resultados/Resultados_Ejercicio3/'+input('Nombre del archivo de salida: ')

config = {
    "approach": approach,
    "numero": num,
    "output_path": output_file
}

# Escribe la configuración en un archivo JSON
with open(path, 'w') as config_file:
    json.dump(config, config_file, indent=4)

print(f"Archivo de configuración '{path}' creado con éxito.")
