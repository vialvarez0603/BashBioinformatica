
#!/usr/bin/env python3

from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
import sys
import json
import logging

logging.basicConfig(filename='errores.log', level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s', filemode='a')
logger = logging.getLogger(__name__)

if len(sys.argv) != 3:
    print("Uso: python script.py <config_file> <mutation_file>")
    sys.exit(1)

archivo_json = sys.argv[1]
path_mutation = sys.argv[2]

if not archivo_json.endswith(".json"):
    logger.error(f"El archivo de entrada '{archivo_json}' no tiene la extensión .json")
    sys.exit(1)

if not path_mutation.endswith(".txt"):
    logger.error(f"El archivo de entrada '{path_mutation}' no tiene la extensión .fasta")
    sys.exit(1)

def extract_seq(archivo):
    with open(archivo, 'r') as file:
        lines = file.readlines()

    secuencia = ""
    for line in lines:
        if not line.startswith(">"):
            secuencia += line.strip()

    return secuencia

def load_config(config_file):
    with open(config_file, 'r') as file:
        config = json.load(file)
    return config

def calculate_melting_temp(primer):
    return mt.Tm_NN(primer)

def is_valid_primer(primer, config):
    gc_content = gc_fraction(primer)
    melting_temp = calculate_melting_temp(primer)
    if config['min_gc_content'] <= gc_content <= config['max_gc_content'] and melting_temp <= config['max_melting_temp']:
        return True
    else:
      return False

def design_primers(transcript, config, number = 5):
    primers = []
    transcript_length = len(transcript)

    for i in range(transcript_length - config['min_length'] + 1):
        for length in range(config['min_length'], config['max_length'] + 1):
            if i + length > transcript_length:
                break
            primer = transcript[i:i + length]
            if is_valid_primer(primer, config) and primer[0] not in "GC" and primer[-1] not in "GC":
                primers.append(primer)
                if len(primers) == number:
                    return primers
    return primers

mutation = extract_seq(path_mutation)
config = load_config(archivo_json)
transcript = mutation
primers = design_primers(transcript, config)

output_file = config['output_path']

if primers:
    print("Primers diseñados:")
    for i,primer in enumerate(primers):
        print(f"Primer {i}:\n\tSecuencia: {primer}, \n\tTm: {calculate_melting_temp(primer):.2f}°C, GC: {gc_fraction(primer) * 100:.2f}%")
else:
    print("No se encontraron primers que cumplan con los criterios.")

with open(output_file, 'w') as f:
    if primers:
        f.write("Primers diseñados:\n")
        for i, primer in enumerate(primers):
            f.write(f"Primer {i}:\n")
            f.write(f"\tSecuencia: {primer},\n")
            f.write(f"\tTm: {calculate_melting_temp(primer):.2f}°C,\n")
            f.write(f"\tGC: {gc_fraction(primer) * 100:.2f}%\n")
    else:
        f.write("No se encontraron primers que cumplan con los criterios.\n")

print(f"Resultados guardados en {output_file}")