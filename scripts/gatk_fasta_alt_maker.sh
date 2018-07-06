#!/bin/sh
fasta="/oasis/projects/nsf/ddp195/dantakli/sperm/somatic_sim/human_g1k_v37.22.fasta"
java -Xmx16g -jar /home/dantakli/bin/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R $fasta -o /oasis/projects/nsf/ddp195/dantakli/sperm/somatic_sim/human_g1k_v37.22.ALT.fasta -L 22 -V $1
