For each gene, we chose the longest sequence from 1 of the 12 Conidae tran- scriptomes as the target sequence. For 421 of these genes, we used the entire length of the sequence as the target sequence, while for the remaining genes, we sliced the target sequences into smaller components based on exon/intron boundaries inferred with EXONERATE v2.2.0 (Slater and Birney 2005) using the Lottia gigantea genome as our reference. EXONERATE v2.2.0 was run under default parameters and under the est2genome model. We chose to use the L. gigantea genome as our reference because it is highly contiguous (scaf- fold N501⁄4 1.87 Mb) and well annotated (Simakov et al. 2013), as a Conidae genome of comparable quality was not available at the time of the bait design. Often, exon/intron boundaries are conserved across fairly divergent taxa and can be used to define exon/intron boundaries (Bi et al. 2012);



https://metazoa.ensembl.org/Lottia_gigantea/Info/Annotation/

la prueba que me gustaría realizar es la siguiente: 
```bash
exonerate --model est2genome $QUERY target.fasta
QUERY=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/californicus/07.Reference/Phuong_etal_2018/01.Assembly/ONLY_Californiconus/Rnaspades_out/transcripts.fasta
```


OR 
https://github.com/yhadevol/Assexon/blob/master/pipeline_scripts/assemble/exonerate_best.pl
```bash
exonerate -m protein2dna queryp/$queryp $contigdir/$subject_sample/$target --score 0 --showvulgar FALSE --showalignment FALSE --subopt FALSE --ryo "sugar: %S %ps\npt: {%Pqs_%Pts_%Ps\t}\n" > $exonerateout/$subject_sample/$out;
```

# Downloading target ref
* 'dna' - unmasked genomic DNA sequences.

* 'dna_rm' - masked genomic DNA.  Interspersed repeats and low complexity regions are detected with the RepeatMasker tool and masked by replacing repeats with 'N's.

* 'dna_sm' - soft-masked genomic DNA. All repeats and low complexity regions have been replaced with lowercased versions of their nucleic base

# Preparing ref

fasta2esd --fasta $REFERENCE --softmask yes --alphabet DNA --output $ESDREF


# Modeling query

-m --model <alignment model>
Specify alignment model type
Supported types:
    ungapped ungapped:trans
    affine:global affine:bestfit affine:local affine:overlap
    est2genome ner protein2dna protein2genome
    protein2dna:bestfit protein2genome:bestfit
    coding2coding coding2genome cdna2genome genome2genome