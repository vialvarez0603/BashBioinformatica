#!/usr/bin/env python3

import json
import sys

# Verifica si se proporcionó un argumento para path
if len(sys.argv) < 2:
    print("Error: Debe proporcionar el nombre del archivo JSON como argumento.")
    sys.exit(1)

# Obtiene el argumento para path
path = sys.argv[1]

print('Determinación de los parámetros de configuración del primer: \nDefault: ')
print('\t Tamaño mínimo del primer: 18p/b')
print('\t Tamaño máximo del primer: 24p/b')
print('\t GC mínimo: 50%')
print('\t GC máximo: 60%')
print('\t Máxima temperatura de melting: 67°C')
print('\t Nombre del archivo de salida: "primers.txt"')

choice = int(input('Para usar la configuración default presionar "0"; Si desea ajustar los parámetros, presione "1":\n'))
while choice != 0 and choice != 1:
  choice = int(input('Vuelva a intentar. Ingrese "0" o "1": '))

if choice == 0:
  mini = 18
  maxi = 24
  min_gc = 50
  max_gc = 60
  max_MT = 67
  output_path = "Resultados/Resultados_Ejercicio5/primers.txt"

else:
  # Define los parámetros de configuración
  mini = int(input('Tamaño mínimo del primer (p/b): '))

  maxi = int(input('Tamaño máximo del primer (p/b): '))
  while maxi < mini:
    maxi = int(input('Error. Debe ser un valor mayor que el previamente ingresado. Vuelva a intentar: '))

  min_gc = float(input('GC mínimo (%): '))
  while min_gc < 0 or min_gc > 100:
    min_gc = float(input('Error. Debe ser un valor entre 0% y 100%. Vuelva a intentar: '))

  max_gc = float(input('GC máximo (%): '))
  while max_gc < min_gc:
    max_gc = float(input('Error. Debe ser un valor mayor que el previamente ingresado. Vuelva a intentar: '))
  while max_gc < 0 or max_gc > 100:
    max_gc = float(input('Error. Debe ser un valor entre 0% y 100%. Vuelva a intentar: '))

  max_MT = float(input('Máxima temperatura de melting (°C): '))

  output_path = 'Resultados/Resultados_Ejercicio5/' + input('Nombre del archivo de salida: ')

config = {
    "min_length": mini,
    "max_length": maxi,
    "min_gc_content": min_gc/100,
    "max_gc_content": max_gc/100,
    "max_melting_temp": max_MT,
    "output_path": output_path
}

# Escribe la configuración en un archivo JSON
with open(path, 'w') as config_file:
    json.dump(config, config_file, indent=4)

print(f"Archivo de configuración {path} creado con éxito.")
