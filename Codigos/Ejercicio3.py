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

def read_alignment(fasta_path):
    alignment = AlignIO.read(fasta_path,"fasta")

    return alignment

def consensus(alignment):
    consensus_seq = []
    num_positions  = alignment.get_alignment_length()

    for column in range(num_positions):
        values = alignment[:,column]
        counts = Counter(values)

        most_common_residue = counts.most_common(1)[0][0]

        if most_common_residue=='-':
            most_common_residue = counts.most_common(2)[1][0]

        consensus_seq.append(most_common_residue)

    return Seq("".join(consensus_seq))

def score(fasta_path, approach, sust = None):
    alignment = read_alignment(fasta_path)
    consensus_seq = consensus(alignment)
    num_positions = alignment.get_alignment_length()
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

def msa2(blast_file, n_blast,output_file):
    sequences = []
    sequence_regex = re.compile(r"<Hsp_hseq>(.*?)</Hsp_hseq>")

    with open(blast_file,"r") as handle:
        #blast_records = NCBIXML.read(result_handle)

        hits_count = 0
        
        for line in handle:
            sequence_match = sequence_regex.search(line)
            if sequence_match:
                sequence = sequence_match.group(1)
                sequences.append(SeqRecord(Seq(sequence), id= f"Hit{hits_count+1}", description = ""))
                hits_count +=1
                if hits_count >= n_blast:
                    break

    temp_fasta = output_path
    
    with open(temp_fasta, "w") as temp_handle:
        for i, sequence in enumerate(sequences):
            temp_handle.write(f">{sequence.id}\n{sequence.seq}\n")

    muscle_cmd = ["muscle", "-in", temp_fasta, "-out", output_path]
    
    try:
        process = subprocess.run(muscle_cmd, check=True, capture_output=True, text=True)
        stdout = process.stdout
        stderr = process.stderr

    except subprocess.CalledProcessError as e:
        logger.error(f"Error ejecutando MUSCLE: {e}")
        logger.error(f"Salida estándar: {e.stdout}")
        logger.error(f"Salida de error: {e.stderr}")
        
        sys.exit(1)
    return

def add_results(fasta_path, score, consensus_seq):
    with open(fasta_path, "r") as original_file:
        original_content = original_file.read()

    new_content = f">Score: {score}\n>Secuencia consenso: {consensus_seq}\n\n{original_content}"
    with open(fasta_path, "w") as modified_file:
        modified_file.write(new_content)
    return

try:
    msa2(resultado_blast, num,output_path)
    score, consensus_seq = score(output_path, app)
    print(f"Score del alineamiento: {score} \nSecuencia consenso: {consensus_seq}")

    add_results(output_path, score, consensus_seq)
except Exception as e:
    logger.error(f"Error en msa2: {e}")
    sys.exit(1)

