
# rutas
EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/cdhit/
export PATH=$PATH:$EXPORT

# correr cdhit con un 99% de similitud (dafault es 90%
# Define CD-HIT-EST parameters
identity=0.95  # 95% identity threshold
word_length=10  # CD-HIT word length parameter
threads=24  # Number of threads for parallel processing

# Run CD-HIT-EST
cd-hit-est -i "${transcriptome}" -o "${OUTFILE}" -c "${identity}" -n "${word_length}" -T "${threads}" -M 20000

echo "CD-HIT-EST completed. Output saved to ${OUTFILE}")

# OR 