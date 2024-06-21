#!/usr/bin/env python3

from Bio.Blast import NCBIWWW #Parametros
import json
import sys
import logging

logging.basicConfig(filename='errores.log', level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s', filemode='a')
logger = logging.getLogger(__name__)

if len(sys.argv) != 3:
    print("Uso: python script.py <config_file> <fasta_interes>")
    sys.exit(1)

config_file = sys.argv[1]
fasta_interes = sys.argv[2]

if not config_file.endswith(".json"):
    logger.error(f"El archivo de entrada '{config_file}' no tiene la extensión .json")
    sys.exit(1)

if not fasta_interes.endswith(".fasta"):
    logger.error(f"El archivo de entrada '{fasta_interes}' no tiene la extensión .fasta")
    sys.exit(1)

try:
    with open(config_file) as f:
        config = json.load(f)
except FileNotFoundError as e:
    logger.error(f"Error al cargar el archivo JSON '{config_file}': {e}")
    sys.exit(1)

tipo = config.get("tipo")
resultados = config.get("resultados")
matriz = config.get("matriz")
e_thres = config.get("e_threshold")
output_path = config.get("output_path")

try:
    with open(fasta_interes, 'r') as fasta_file:
        fasta_data = fasta_file.read().strip()
        print(f"fasta_data: \n{fasta_data}")
except FileNotFoundError as e:
    logger.error(f"Error al cargar el archivo FASTA '{fasta_interes}': {e}")
    sys.exit(1)

try:
    result = NCBIWWW.qblast('blastp','swissprot',fasta_data,descriptions=resultados,alignments=resultados,hitlist_size=resultados,matrix_name=matriz,expect=e_thres)
    blast=open(output_path,'w')
    blast.write(result.read())
    blast.close()
except Exception as e:
    logger.error(f"Error ejecutando el BLAST: {e}")
    sys.exit(1)

