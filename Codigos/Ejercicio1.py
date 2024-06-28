
#!/usr/bin/env python3

import Bio
from Bio import SeqIO, AlignIO, motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from collections import Counter
import numpy as np
import json
import logging
import sys

def extract_orfs(record):
    orf_records = []
    orf_lengths = []
    for frame in range(3):
        seq = record.seq[frame:]
        protein_seq = seq.translate()
        orf_record = SeqRecord(protein_seq, id=f"{record.id}_frame{frame+1}", description=f"ORF_{frame+1}")
        orf_records.append(orf_record)
        orf_lengths.append(len(orf_record))

    return orf_records, orf_lengths

def main(input_file, output_file):
    records = SeqIO.parse(input_file, "genbank")

    with open(output_file, "w") as output_handle:
        for record in records:
            complemento = record.seq[::-1]
            complemento_record = SeqRecord(complemento, id=f"{record.id}_inverse", description="Inverso")

            orf_records = extract_orfs(record)[0] + extract_orfs(complemento_record)[0]
            SeqIO.write(orf_records, output_handle, "fasta")

def calc_protein_lengths(ORF_seq):
    lengths = []
    protein_length = 0
    for aa in ORF_seq:
        if aa != "*":
            protein_length += 1
        else:
            lengths.append(protein_length)
            protein_length = 0

    return lengths

def protein_lengths_per_ORF(fasta_file):
    orf_proteins = []
    with open(fasta_file) as f:
        orf_seq = ""
        orf_id = None
        for line in f:
            if line.startswith(">"):
                if orf_id is not None:
                    lengths = calc_protein_lengths(orf_seq)
                    lengths = sorted(lengths, reverse=True)
                    orf_proteins.append((orf_id, lengths))
                orf_id = line.strip()
                orf_seq = ""
            else:
                orf_seq += line.strip()

        if orf_id is not None:
            lengths = calc_protein_lengths(orf_seq)
            lengths = sorted(lengths, reverse=True)
            orf_proteins.append((orf_id, lengths))

    return orf_proteins

def choose_orf(orfs, object_length):
    min_difs = []
    for orf_id, lengths in orfs:
        min_dif = np.min(np.abs(np.ones(len(lengths))*object_length-lengths))
        min_difs.append([orf_id,min_dif])

    lista_ordenada = sorted(min_difs, key=lambda x: x[1])

    return lista_ordenada

def buscar_fasta(ID,file):
    with open(file) as f:
        orf_id = None
        fasta_interes=""
        for line in f:
            if line.startswith(">"):
                orf_id = line.strip()
            elif orf_id == ID:
                fasta_interes += line.strip()
    return fasta_interes

logging.basicConfig(filename='errores.log', level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_from_json(config_file):
    try:
        with open(config_file) as f:
            config = json.load(f)
    except FileNotFoundError as e:
        logger.error(f"Error al cargar el archivo JSON '{config_file}': {e}")
        sys.exit(1)

    input_file = config.get("input_file")
    output_file_fas = config.get("output_file")
    file_of_interest = config.get("file_of_interest")
    longitud = config.get("longitud")

    # Verificación de extensiones
    if not input_file.endswith(".gb"):
        logger.error(f"El archivo de entrada '{input_file}' no tiene la extensión .gb")
        sys.exit(1)
    if not output_file_fas.endswith(".fasta"):
        logger.error(f"El archivo de salida '{output_file_fas}' no tiene la extensión .fasta")
        sys.exit(1)
    if not file_of_interest.endswith(".fasta"):
        logger.error(f"El archivo de interés '{file_of_interest}' no tiene la extensión .fasta")
        sys.exit(1)

    main(input_file, output_file_fas)
    return output_file_fas, file_of_interest, longitud

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python script.py <config_file>")
        sys.exit(1)

    config_file = sys.argv[1]
    fileof_interest = sys.argv[2]
    output_file_fas, file_of_interest, length = run_from_json(config_file)

    longitud_proteinas = protein_lengths_per_ORF(output_file_fas)
    for orf_id, lengths in longitud_proteinas:
        print(f"ORF: {orf_id}")
        print(f"Longitudes de las proteínas codificadas: {lengths}")


    orf_ordenada = choose_orf(longitud_proteinas,length)

    mini = np.infty
    mini_orf = ''
    for orf, dif in orf_ordenada:
        print(f"ORF: {orf} - Diferencia mínima: {dif} AA")
        if dif < mini:
            mini = dif
            mini_orf = orf
    
    print(f"ORF seleccionado: {mini_orf}")
    fasta_final = buscar_fasta(mini_orf,output_file_fas)
    with open(file_of_interest, "w") as f:
        f.write(fasta_final + "\n")
