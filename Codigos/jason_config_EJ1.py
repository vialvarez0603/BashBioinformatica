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
print('\tNombre del archivo Genbank: "input_genbank.gb"')
print('\tNombre del archivo FASTA resultante: "output_FASTA.fasta"')
print('\tLongitud de la proteína de interés: 140aa')
print('\tNombre del archivo FASTA correspondiente al ORF seleccionado: "ORF_FASTA.fasta"')

choice = int(input('Para usar la configuración default presionar "0"; Si desea ajustar los parámetros, presione "1":\n'))
while choice != 0 and choice != 1:
    choice = int(input('Vuelva a intentar. Ingrese "0" o "1": '))

if choice == 0:
    input_file = "input_genbank.gb"
    output_file_fas = "Resultados/Resultados_Ejercicio1/output_FASTA.fasta"
    longitud = 140
    archivo_interes = "Resultados/Resultados_Ejercicio1/ORF_FASTA.fasta"
else:
    # Define los parámetros de configuración
    input_file = input('Nombre del archivo Genbank: ')
    output_file_fas = 'Resultados/Resultados_Ejercicio1/' + input('Nombre del archivo FASTA resultante: ')
    longitud = int(input('Longitud de la proteína de interés: '))
    while longitud <= 0:
        longitud = int(input('Error. Debe ser un valor mayor a 0. Vuelva a intentar: '))
    archivo_interes = 'Resultados/Resultados_Ejercicio1/' + input('Nombre del archivo FASTA correspondiente al ORF seleccionado: ')

config = {
    "input_file": input_file,
    "output_file": output_file_fas,
    "longitud": longitud,
    "file_of_interest": archivo_interes
}

# Escribe la configuración en un archivo JSON
with open(path, 'w') as config_file:
    json.dump(config, config_file, indent=4)

print(f"Archivo de configuración '{path}' creado con éxito.")
