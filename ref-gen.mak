PREFIX=hs38
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
PATH_TOOLS=.
N_THREADS=8

.PHONY:all genome bwa_index bowtie2_index tgz

all:genome bwa_index bowtie2_index

genome:$(PREFIX)/$(PREFIX).fa $(PREFIX)/$(PREFIX).fa.fai $(PREFIX)/$(PREFIX).dict

bwa_index:$(PREFIX)/$(PREFIX).fa.bwt

bowtie2_index:$(PREFIX)/$(PREFIX).fa.1.bt2

tgz:$(PREFIX).genome.tgz $(PREFIX).bwa_index.tgz $(PREFIX).bowtie2_index.tgz

$(PREFIX).mkdir:
		mkdir -p $(PREFIX) && touch $@

$(PREFIX)/$(PREFIX).fa:$(PREFIX).mkdir
		curl -L --retry 10 $(URL) | $(PATH_TOOLS)/seqtk seq -Ul80 > $@

$(PREFIX)/$(PREFIX).fa.fai:$(PREFIX)/$(PREFIX).fa
		$(PATH_TOOLS)/samtools faidx $<

$(PREFIX)/$(PREFIX).dict:$(PREFIX)/$(PREFIX).fa
		java -jar $(PATH_TOOLS)/picard.jar CreateSequenceDictionary -R $< -O $@ --URI $(URL)

$(PREFIX)/$(PREFIX).fa.bwt:$(PREFIX)/$(PREFIX).fa
		$(PATH_TOOLS)/bwa index $<

$(PREFIX)/$(PREFIX).fa.1.bt2:$(PREFIX)/$(PREFIX).fa
		$(PATH_TOOLS)/bowtie2-build-s $< $<

$(PREFIX).genome.tgz:$(PREFIX)/$(PREFIX).fa $(PREFIX)/$(PREFIX).fa.fai $(PREFIX)/$(PREFIX).dict
		tar -cf - $^ | $(PATH_TOOLS)/pigz -c -p$(N_THREADS) > $@

$(PREFIX).bwa_index.tgz:$(PREFIX)/$(PREFIX).fa.bwt
		tar -cf - $(PREFIX)/$(PREFIX).fa.amb $(PREFIX)/$(PREFIX).fa.ann $(PREFIX)/$(PREFIX).fa.bwt $(PREFIX)/$(PREFIX).fa.pac $(PREFIX)/$(PREFIX).fa.sa | $(PATH_TOOLS)/pigz -c -p$(N_THREADS) > $@

$(PREFIX).bowtie2_index.tgz:$(PREFIX)/$(PREFIX).fa.1.bt2
		tar -cf - $(PREFIX)/$(PREFIX).fa.?.bt2 $(PREFIX)/$(PREFIX).fa.rev.?.bt2 | $(PATH_TOOLS)/pigz -c -p$(N_THREADS) > $@
