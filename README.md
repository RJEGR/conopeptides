# De novo transcriptome assembly for the specie Californiconus californicus
# Author
Ricardo Gomez-Reyes
# Description
<div align="justify">
The justified description of the proyect will be here
</div>

Here good introduction to methods for conopeptides [identification](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-60#Sec15)

Continue in 
/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/californicus/02.Assembly

Here we propose an alternative representation for each gene, which we refer to as a superTranscript. SuperTranscripts contain the sequence of all exons of a gene without redundancy (Figure 1A). They can be constructed from any set of transcripts including de novo assemblies and we have developed a python program to build them called Lace (available from https://github.com/Oshlack/Lace/wiki). Lace works by building a splice graph[13] for each gene, then topologically sorting the graph using Kahn’s algorithm[14] (Figure 1B). Building superTranscripts is a simple post-assembly step that promises to unlock numerous analytical approaches for non-model organisms.

A powerful, new application of superTranscripts is merging transcriptomes from a variety of sources. 