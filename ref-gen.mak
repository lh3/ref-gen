PREFIX=hs38
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
PATH_TOOLS=.

.PHONY:all

all:$(PREFIX)/$(PREFIX).fa $(PREFIX)/$(PREFIX).fa.fai $(PREFIX)/$(PREFIX).dict $(PREFIX)/$(PREFIX).fa.bwt $(PREFIX)/$(PREFIX).fa.1.bt2

$(PREFIX).mkdir:
		mkdir -p $(PREFIX) && touch $@

$(PREFIX)/$(PREFIX).fa:$(PREFIX).mkdir
		curl -L $(URL) | gzip -dc | $(PATH_TOOLS)/seqtk seq -Ul80 > $@ || echo done

$(PREFIX)/$(PREFIX).fa.fai:$(PREFIX)/$(PREFIX).fa
		$(PATH_TOOLS)/samtools faidx $<

$(PREFIX)/$(PREFIX).dict:$(PREFIX)/$(PREFIX).fa
		java -jar $(PATH_TOOLS)/picard.jar CreateSequenceDictionary -R $< -O $@ --URI $(URL)

$(PREFIX)/$(PREFIX).fa.bwt:$(PREFIX)/$(PREFIX).fa
		$(PATH_TOOLS)/bwa index $<

$(PREFIX)/$(PREFIX).fa.1.bt2:$(PREFIX)/$(PREFIX).fa
		$(PATH_TOOLS)/bowtie2-build-s $< $<
