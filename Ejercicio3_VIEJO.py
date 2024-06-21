#!/usr/bin/env python3

import Bio
from Bio import SeqIO, AlignIO, motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from Bio.Blast import NCBIXML
from io import StringIO
from collections import Counter
import numpy as np
import re
import os
import subprocess
import json
import sys
import logging

logging.basicConfig(filename='errores.log', level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s', filemode='a')
logger = logging.getLogger(__name__)

if len(sys.argv) != 3:
    print("Uso: python script.py <config_file> <blast_file>")
    sys.exit(1)

config_file = sys.argv[1]
resultado_blast = sys.argv[2]

if not config_file.endswith(".json"):
    logger.error(f"El archivo de entrada '{config_file}' no tiene la extensión .json")
    sys.exit(1)

if not resultado_blast.endswith(".txt"):
    logger.error(f"El archivo de entrada '{resultado_blast}' no tiene la extensión .txt")
    sys.exit(1)

try:
    with open(config_file) as f:
        config = json.load(f)
except FileNotFoundError as e:
    logger.error(f"Error al cargar el archivo JSON '{config_file}': {e}")
    sys.exit(1)

app = config.get("approach")
num = config.get("numero")
output_path = config.get("output_path")

def consensus(aligned_motif):
    consensus_seq = []
    alignment = aligned_motif.alignment
    for column in range(len(alignment[0])):
      values = [record[column] for record in alignment]
      counts = Counter(values)

      most_common_residue = counts.most_common(1)[0][0]

      if most_common_residue=='-':
        most_common_residue = counts.most_common(2)[1][0]

      consensus_seq.append(most_common_residue)

    return Seq("".join(consensus_seq))

def score(aligned_motif, approach, sust = None):
    alignment = aligned_motif.alignment
    consensus_seq = consensus(aligned_motif)
    num_positions = len(alignment[0])
    num_secuencias = len(alignment)

    if sust == None:
      sust_mat = [1 if i == j else 0 for j in range(num_positions) for i in range(num_secuencias)]

    if approach == 1:
      scores = []
      for column in range(num_positions):
        count = num_secuencias - sum(1 for row in alignment if row[column] == consensus_seq[column])
        scores.append(count)

      score = -sum(scores)

    elif approach == 2:
      scores = []
      for column in range(num_positions):
        count = sum(1 for row in alignment if row[column] == consensus_seq[column])
        tot = count*2 - (num_secuencias-count)
        scores.append(tot)

      score = sum(scores)

    else:
      score = None

    return score, consensus_seq

def msa2(blast_file, n_blast, approach = 1):
    sequences = []

    with open(blast_file,"r") as result_handle:
        #blast_records = NCBIXML.read(result_handle)

        hits_count = 0
        '''for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                if hits_count >= n_blast:
                    break
                sequences.append(SeqRecord(Seq(hsp.sbjct), id=f"Hit{hits_count+1}", description=""))
                hits_count += 1
            if hits_count >= n_blast:
                break'''
        for line in handle:
            sequence_match = sequence_regex.search(line)
            if sequence_match:
                sequence = sequence_match.group(1)
                sequences.append(SeqRecord(Seq(sequence)), id= f"Hit{hits_count+1}", description = ""))
                hits_count +=1
                if hits_count >= n_blast:
                    break

    temp_fasta = "temp.fasta"
    temp_output = "out.fasta"
    with open(temp_fasta, "w") as temp_handle:
        for i, sequence in enumerate(sequences):
            temp_handle.write(f">{sequence.id}\n{sequence.seq}\n")
       #SeqIO.write(sequences, temp_handle, "fasta")

    muscle_cmd = ["muscle", "-in", temp_fasta, "-out", temp_output]
    try:
        process = subprocess.run(muscle_cmd, check=True, capture_output=True, text=True)
        stdout = process.stdout
        stderr = process.stderr

    except subprocess.CalledProcessError as e:
        logger.error(f"Error ejecutando MUSCLE: {e}")
        logger.error(f"Salida estándar: {e.stdout}")
        logger.error(f"Salida de error: {e.stderr}")
        
        sys.exit(1)
    
    '''with open(temp_output, "r") as f:
        print("Contenido de out.fasta")
        print(f.read())'''

    try:
        alignment = AlignIO.read(temp_output, "fasta").alignmente
        motif = motifs.Motif(alignment=alignment)
        alignment_score, consensus_seq = score(motif, approach)

        os.remove(temp_fasta)
        os.remove(temp_output)

        return motif, alignment, consensus_seq, alignment_score

    except Exception as e:
        logger.error(f"Error en la lectura de alineamiento: {e}")
        sys.exit(1)
try:
    motif, alineamiento, consenso, score_tot = msa2(resultado_blast, num, approach=app)
except Exception as e:
    logger.error(f"Error en msa2: {e}")
    sys.exit(1)

print(f"Alineamiento: \n{alineamiento}")

print(f"Secuencia consenso: \n{consenso}\n")
print(f"Score del alineamiento: \n{score_tot}\n")

try:
    with open(output_path, 'w') as f:
        f.write(f"Alineamiento: \n{alineamiento}\n\n")
        f.write(f"Secuencia consenso: \n{consenso}\n\n")
        f.write(f"Score del alineamiento: \n{score_tot}\n")
    print(f"Resultados guardados en {output_path}")
except IOError as e:
    logger.error(f"Error al escribir en el archivo '{output_path}': {e}")
    sys.exit(1)
