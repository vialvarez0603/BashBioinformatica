log_file="errores.log"

if [ -f "errores.log" ]; then
	rm "errores.log"
fi

chmod +x Codigos/muscle5.1.linux_intel64

#EJERCICIO 1
#python3 Codigos/jason_config_EJ1.py Resultados/Resultados_Ejercicio1/config_EJ1.json

python3 Codigos/Ejercicio1.py Resultados/Resultados_Ejercicio1/config_EJ1.json Resultados/Resultados_Ejercicio1/ORF_FASTA.fasta

#EJERCICIO 2
#python3 Codigos/jason_config_EJ2.py Resultados/Resultados_Ejercicio2/config_EJ2.json

python3 Codigos/Ejercicio2.py Resultados/Resultados_Ejercicio2/config_EJ2.json Resultados/Resultados_Ejercicio1/ORF_FASTA.fasta

#EJERCICIO 3
#python3 Codigos/jason_config_EJ3.py Resultados/Resultados_Ejercicio3/config_EJ3.json

python3 Codigos/Ejercicio3.py Resultados/Resultados_Ejercicio3/config_EJ3.json Resultados/Resultados_Ejercicio2/RESULTADO_BLAST.txt

#EJERCICIO 4
# Definir los nombres de los archivos de entrada y salida

archivo_fasta="gen.fasta"
salida_fasta="orf_output.fasta"
salida_aa="aa_output.fasta"
archivo_orf_fasta="Resultados/Resultados_Ejercicio1/ORF_FASTA.fasta"

# Solicita el archivo del gen en formato FASTA
while true; do
    echo -e "\033[1;37mVerificando la existencia del archivo del gen, en formato fasta:\033[0m" 

    # Verifica si el archivo existe
    if [ ! -f "$archivo_fasta" ]; then
        echo -e "\033[1;31mEl archivo no existe. El archivo que debe ingresar se llama gen.fasta, por favor inténtelo de nuevo.\033[0m"|tee -a "$log_file"
        exit 1  # Sale del script si el archivo no existe
    else
        break  
    fi
done

# Define la ruta base para los resultados
nombre_base="."

# Crea la carpeta para los resultados del Ejercicio4. Si ya existe, se borra y se vuelve a crear
resultados_dir="$nombre_base/Resultados/Resultados_Ejercicio4"
if [ -d "$resultados_dir" ]; then
  rm -rf "$resultados_dir"
fi
mkdir -p "$resultados_dir"

# Generar los ORFs y guardarlos en el archivo de salida predefinido
getorf -sequence "$archivo_fasta" -outseq "$resultados_dir/$salida_fasta"

# Traducir las secuencias de nucleótidos a aminoácidos y guardarlas en el archivo de salida predefinido
transeq -sequence "$archivo_fasta" -outseq "$resultados_dir/$salida_aa"

# Verificar la existencia del archivo de ORF
while true; do
    echo -e "\033[1;37mVerificando la existencia del archivo de ORF fasta:\033[0m"

    # Verifica si el archivo existe
    if [ ! -f "$archivo_orf_fasta" ]; then
        echo -e "\033[1;31mEl archivo no existe. El archivo que debe ingresar se llama ORF_correcto.fasta, por favor inténtelo de nuevo.\033[0m" | tee -a "$log_file"
        exit 1  # Sale del script si el archivo no existe
    else
        break  
    fi
done

# Realiza el análisis de motivos usando la base de datos PROSITE y guarda el resultado en un archivo de salida predefinido
#/bin/prosextract -sequence "$archivo_orf_fasta" -outfile "$resultados_dir/motivos.fasta"

prosite_dir=$(pwd)
prositedir="$prosite_dir/prositedir"
prosite_data="$prositedir/PROSITE"

mkdir -p "$prosite_data"

if [ -f "$prosite_dir/prosite.dat" ] && [ -f "$prosite_dir/prosite.doc" ]; then
	cp "$prosite_dir/prosite.dat" "$prosite_data/"
	cp "$prosite_dir/prosite.doc" "$prosite_data/"

	chmod 777 "$prosite_data/prosite.dat"
	chmod 777 "$prosite_data/prosite.doc"
else
	echo -e "\033[1;31mNo se encontraron prosite.dat o prosite.doc en el directorio actual. Aeg[urese de que ambos archivos est[en presentes.\033[0m" | tee -a "$log_file"
	exit 1
fi

if [ -f "$prosite_data/prosite.dat" ] && [ -f "$prosite_data/prosite.doc" ]; then
	echo -e "\033[1;32mArchivos prosite.dat y prosite.doc copiados exitosamente a $prosite_data.\033[0m"
else
	echo -e "\033[1;32mFallo al copiar prosite.dat o prosite.doc a $prosite_data.\033[0m"
	exit 1
fi

export EMBOSS_DATA="$prositedir"

prosextract -prositedir "$prositedir" 2>>"$log_file"
patmatmotifs -sequence "$archivo_orf_fasta" -outfile "$resultados_dir/motivos.txt" 2>>"$log_file"
#patmatmotifs -outfile "$resultados_dir/motivos.txt" "$archivo_orf_fasta" 2>>"$log_file"

echo "Análisis completado. Resultados guardados en $resultados_dir"

#EJERCICIO 5
#python3 Codigos/jason_config_EJ5.py Resultados/Resultados_Ejercicio5/config_EJ5.json

python3 Codigos/Ejercicio5.py Resultados/Resultados_Ejercicio5/config_EJ5.json mutacion.txt
	

