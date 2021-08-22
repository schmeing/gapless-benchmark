# gapless-benchmark
Snakemake file to reproduce the benchmark in the gapless paper

The datasets, references and software has to be downloaded from their original source.
The E. coli data can be downloaded [Public Health England reference collections](https://www.sanger.ac.uk/resources/downloads/bacteria/nctc/) and the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/SRR3191692), the dolphin data from the [Vertebrate Genomes Project](https://vgp.github.io/genomeark/Tursiops_truncatus/) and the human data from the [T2T consortium](https://github.com/marbl/CHM13)

The files not included, but necessary to reproduce the figures in the paper are:
```
input/truncatus/reference/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna
input/truncatus/PacBio_CLS/ (full folder)
input/truncatus/chromium/ (full folder)
input/ecoli/SRR3191692/SRR3191692_1.fq.gz
input/ecoli/SRR3191692/SRR3191692_2.fq.gz
input/ecoli/reference/ERS764956.gff
input/ecoli/ERR1036235/ (full folder)
input/human/T2T_10X_NovaSeq/ (full folder)
input/human/PRJNA559484_Nanopore/rel7.fastq.gz
input/human/reference/chm13.draft_v1.1.fasta.gz
input/human/PRJNA530776_SequelII_HiFi/ (full folder)
```
