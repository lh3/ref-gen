## Introduction

This repo provides the URLs of several human reference genomes. These genomes
were adapted from the original GenBank verions and work better for read
alignment and variant calling. All of them include [unlocalized][GRC-def] and
unplaced contigs, use [rCRS][rCRS] for mitochondrion and have pseudoautosomal
regions (PARs) on Y chromosome hard masked.

To download these genomes, see links in [URL.txt](URL.txt). The genomes are
also "permanently" avaialble at doi:[10.5281/zenodo.8045373][zenodo] along with
their BWA and Bowtie2 indices. Command lines for generating secondary files can
be found in [ref-gen.mak](ref-gen.mak). Note that downloading from Zenodo can
be slow. I will look for a faster host in future.

## Description

The following table gives a brief description of these genomes:

|Name   |Source |Brief description|
|:------|:------|:----------|
|hs37   |GRCh37 |1000 Genomes Project (1000G) reference in 2010|
|hs37d5 |GRCh37 |hs37 with decoy (recommended GRCh37)|
|hs38   |GRCh38 |GRCh38 no-alt analysis set (recommended GRCh38)|
|hs38DH |GRCh38 |GRCh38 with ALT, decoy and HLA genes (not recommended)|
|chm13v2|CHM13v2|T2T CHM13 plus HG002 chrY (recommended CHM13)|

In more detail:

* **hs37**. This was the reference genome used by the 1000 Genomes Project at
  the early stage. It was derived from the Ensembl primary genome but had
  mitochondrial DNA replaced with rCRS.

* **hs37d5**. This was generated on top of **hs37** and used by 1000G in 2011.
  The only difference is the addition of decoy sequences and the EBV genome.
  This version has been found to give better variant calls. Recommended version
  if you use GRCh37.

* **hs38**. This is the so-called GRCh38 no-ALT analysis set. Recommended if
  you use GRCh38.

* **hs38DH**. This genome includes GRCh38 alternate contigs, GRCh38 decoy
  sequences and HLA alleles. As GRCh38 is more complete than GRCh37, GRCh38
  decoy sequences are not as important as GRCh37 decoy. Furthermore, most tools
  would not work well with this version as they are not ALT-aware. Improper
  use of this genome would greatly hurt variant sensitivity. Not recommended
  unless you know what you are doing.

I was involved in discussions on generating **hs37** and **hs38**. I created
the final FASTA files of **hs37d5**, **hs38DH** and **chm13v2**. Of course,
the original genome sequences were produced by others. If you are confused by
**hs38DH**, you may blame me for that.

Note that the description above is consistent with my [old blog post][blog] on
"Which human reference genome to use?" This time I added **chm13v2** as well as
download links to packaged sequences and indices.

[GRC]: https://www.ncbi.nlm.nih.gov/grc
[rCRS]: https://www.mitomap.org/MITOMAP/HumanMitoSeq
[GRC-def]: https://www.ncbi.nlm.nih.gov/grc/help/definitions/
[zenodo]: https://zenodo.org/record/8045374
[hs37d5]: https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/
[blog]: http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
