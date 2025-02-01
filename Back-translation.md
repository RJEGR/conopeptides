# Sequence Translation 



Back-translation is used to predict the possible nucleic acid sequence that a specified peptide sequence has originated from.

EMBOSS Backtranambig [here](https://www.ebi.ac.uk/jdispatcher/st)
+ Purpose: Backtranambig, on the other hand, back-translates a protein sequence into an ambiguous nucleotide sequence. This means it accounts for the degeneracy of the genetic code, allowing for multiple possible nucleotide sequences for each amino acid.
+ Use Case: This tool is beneficial when you want to explore all possible nucleotide sequences that could encode a given protein, especially in cases where codon preferences are not known or when considering evolutionary studies

Use Backtranseq for specific nucleotide sequences when codon usage is known, and use Backtranambig when you want to account for multiple possible sequences due to codon degeneracy.

Due to web app limit sequences to 500, I will to locally install and run.

```bash
EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/EMBOSS-6.6.0/emboss
export PATH=$PATH:$EXPORT

WD=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/DataBases/Conopeptides
cd $WD

# backtranambig -h
backtranambig -sequence conoserver_protein.curated.fa -outfile conoserver_protein.curated_backtranambig.fa

# 
# srun backtranseq -sequence conoserver_protein.curated.fa -outfile conoserver_protein.curated_backtranseq.fa
```

The codon usage table is read by default from "Ehum.cut" in the 'data/CODONS' directory of the EMBOSS distribution. 
/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/EMBOSS-6.6.0/emboss/data/CODONS/
If the name of a codon usage file is specified on the command line, then this file will first be searched for in the current directory and then in the 'data/CODONS' directory of the EMBOSS distribution


Other options 
+ [DNAChisel](https://edinburgh-genome-foundry.github.io/DnaChisel/_modules/dnachisel/biotools/sequences_operations.html#reverse_translate)
+ [codon_tools](https://github.com/jocelynnpearl/codon_tools/tree/master)