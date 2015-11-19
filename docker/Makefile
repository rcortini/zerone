URL= http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncode
TWOBIT= hg19.2bit
FASTA= hg19_unmasked.fasta
IDX= hg19_unmasked.gem
Z= ../zerone/
B= ../BayesPeak/
J= ../JAMM/
M= ../MACS/

all: curves.pdf
	md5sum -c datasets.md5

wgEncodeUwTfbsK562InputStdRawDataRep1.fastq.gz:
	# Control of CTCF and H3K36me3.
	wget $(URL)UwTfbs/$@

wgEncodeUwTfbsK562CtcfStdRawDataRep%.fastq.gz:
	# Targets of CTCF.
	wget $(URL)UwTfbs/$@

wgEncodeUwHistoneK562H3k36me3StdRawDataRep%.fastq.gz:
	# Targets of H3K36me3.
	wget $(URL)UwHistone/$@

wgEncodeSydhTfbsK562InputRawData.fastq.gz:
	# Control of Pol2.
	wget $(URL)SydhTfbs/$@

wgEncodeSydhTfbsK562Pol2RawDataRep%.fastq.gz:
	# Targets of Pol2.
	wget $(URL)SydhTfbs/$@

$(TWOBIT):
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/$@

twoBitToFa:
	wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/$@
	chmod +x $@

$(FASTA): $(TWOBIT) twoBitToFa
	./twoBitToFa $< $@

$(IDX): $(FASTA)
	gem-indexer --complement emulate -i $< -o $(subst .gem,,$@) -T 4

%.map: %.fastq.gz $(IDX)
	zcat $< | \
	gem-mapper -I $(IDX) --unique-mapping -q ignore -m 2 -T 4 \
	-o $(subst .map,,$@)

%_no_random.map: %.map ignore_random_chrom.py
	./ignore_random_chrom.py $< > $@

$BPol%_output.bed: $BBayesPeak2bed.py


$ZCtcf.bed: wgEncodeUwTfbsK562InputStdRawDataRep1_no_random.map wgEncodeUwTfbsK562CtcfStdRawDataRep1_no_random.map wgEncodeUwTfbsK562CtcfStdRawDataRep2_no_random.map
	$Zzerone -l \
	--mock wgEncodeUwTfbsK562InputStdRawDataRep1_no_random.map \
	--chip wgEncodeUwTfbsK562CtcfStdRawDataRep1_no_random.map \
	--chip wgEncodeUwTfbsK562CtcfStdRawDataRep2_no_random.map \
	> $@

$ZH3k36me3.bed: wgEncodeUwTfbsK562InputStdRawDataRep1_no_random.map wgEncodeUwHistoneK562H3k36me3StdRawDataRep1_no_random.map wgEncodeUwHistoneK562H3k36me3StdRawDataRep2_no_random.map
	$Zzerone -l \
	--mock wgEncodeUwTfbsK562InputStdRawDataRep1_no_random.map \
	--chip wgEncodeUwHistoneK562H3k36me3StdRawDataRep1_no_random.map \
	--chip wgEncodeUwHistoneK562H3k36me3StdRawDataRep2_no_random.map \
	> $@

$ZPol2.bed: wgEncodeSydhTfbsK562InputRawData_no_random.map wgEncodeSydhTfbsK562Pol2RawDataRep1_no_random.map wgEncodeSydhTfbsK562Pol2RawDataRep2_no_random.map
	$Zzerone -l \
	--mock wgEncodeSydhTfbsK562InputRawData_no_random.map \
	--chip wgEncodeSydhTfbsK562Pol2RawDataRep1_no_random.map \
	--chip wgEncodeSydhTfbsK562Pol2RawDataRep2_no_random.map \
	> $@

#CTCF.txt:
#H3K36me3.txt:
#Pol2.txt:
curves.pdf: $BPol1_output.bed $BPol2_output.bed $ZCtcf.bed $ZH3k36me3.bed $ZPol2.bed
	R -f plot_curves.R
