URL= http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncode
TWOBIT= hg19.2bit
FASTA= hg19_unmasked.fasta
IDX= hg19_unmasked.gem

B= ../BayesPeak/
J= ../JAMM/
M= ../MACS/
Z= ../zerone/

BayesPeak= $BBayesPeak.R $BBayesPeak2bed.py
CHRSIZE= $Jhg19_chrom_sizes

Uw_input= wgEncodeUwTfbsK562InputStdRawDataRep1
Pol2_input= wgEncodeSydhTfbsK562InputRawData
Ctcf_target= wgEncodeUwTfbsK562CtcfStdRawDataRep
H3k4me3_target= wgEncodeUwHistoneK562H3k4me3StdRawDataRep
H3k36me3_target= wgEncodeUwHistoneK562H3k36me3StdRawDataRep
Pol2_target= wgEncodeSydhTfbsK562Pol2RawDataRep
Ctcf_all= $(Uw_input) $(Ctcf_target)1 $(Ctcf_target)2
H3k4me3_all= $(Uw_input) $(H3k4me3_target)1 $(H3k4me3_target)2
H3k36me3_all= $(Uw_input) $(H3k36me3_target)1 $(H3k36me3_target)2
Pol2_all= $(Pol2_input) $(Pol2_target)1 $(Pol2_target)2

PROF= Ctcf H3k4me3 H3k36me3 Pol2
BOUT= $(addprefix $B,$(addsuffix 1.bed,$(PROF))) \
      $(addprefix $B,$(addsuffix 2.bed,$(PROF)))
JOUT= $(addprefix $J,$(addsuffix .bed,$(PROF)))
MOUT= $(addprefix $M,$(addsuffix 1.bed,$(PROF))) \
      $(addprefix $M,$(addsuffix 2.bed,$(PROF)))
ZOUT= $(addprefix $Z,$(addsuffix .bed,$(PROF)))

.SECONDARY:

all: $(BOUT) $(JOUT) $(MOUT) $(ZOUT)
	md5sum -c datasets.md5

$(Uw_input).fastq.gz:
	wget $(URL)UwTfbs/$@

$(Ctcf_target)%.fastq.gz:
	wget $(URL)UwTfbs/$@

$(H3k4me3_target)%.fastq.gz:
	wget $(URL)UwHistone/$@

$(H3k36me3_target)%.fastq.gz:
	wget $(URL)UwHistone/$@

$(Pol2_input).fastq.gz:
	wget $(URL)SydhTfbs/$@

$(Pol2_target)%.fastq.gz:
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

fetchChromSizes:
	wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/$@

$(CHRSIZE): fetchChromSizes
	bash $< hg19 > $@

%.map: %.fastq.gz $(IDX) ignore_random_chrom.py
	zcat $< | \
	gem-mapper -I $(IDX) --unique-mapping -q ignore -m 2 -T 4 | \
	./ignore_random_chrom.py > $@

%.sam: %.map
	gem-2-sam -I $(IDX) -l -q 'offset-64' -i $< > $@

wgEncode%.bed: wgEncode%.sam
	sam2bed $< $@

$BCtcf%.bed: $(Ctcf_target)%.bed $(Uw_input).bed $(BayesPeak)
	Rscript $BBayesPeak.R $< $(Uw_input).bed
	$BBayesPeak2bed.py tmp > $@
	rm tmp

$BH3k4me3%.bed: $(H3k4me3_target)%.bed $(Uw_input).bed $(BayesPeak)
	Rscript $BBayesPeak.R $< $(Uw_input).bed
	$BBayesPeak2bed.py tmp > $@
	rm tmp

$BH3k36me3%.bed: $(H3k36me3_target)%.bed $(Uw_input).bed $(BayesPeak)
	Rscript $BBayesPeak.R $< $(Uw_input).bed
	$BBayesPeak2bed.py tmp > $@
	rm tmp

$BPol2%.bed: $(Pol2_target)%.bed $(Pol2_input).bed $(BayesPeak)
	Rscript $BBayesPeak.R $< $(Pol2_input).bed
	$BBayesPeak2bed.py tmp > $@
	rm tmp

$JCtcf.bed: $(addsuffix .bed,$(Ctcf_all)) $JJAMM.sh $(CHRSIZE)
	mkdir -p $JCtcf/input $JCtcf/target
	cp $(Uw_input).bed $JCtcf/input/$(Uw_input).bed
	cp $(Ctcf_target)1.bed $JCtcf/target/$(Ctcf_target)1.bed
	cp $(Ctcf_target)2.bed $JCtcf/target/$(Ctcf_target)2.bed
	bash $JJAMM.sh \
	    -g $(CHRSIZE) \
	    -o $JCtcf \
	    -s $JCtcf/target \
	    -c $JCtcf/input

$JH3k4me3.bed: $(addsufifix .bed,$(H3k4me3_all)) $JJAMM.sh $(CHRSIZE)
	mkdir -p $JH3k4me3/input $JH3k4me3/target
	cp $(Uw_input).bed $JH3k4me3/input/$(Uw_input).bed
	cp $(H3k4me3_target)1.bed $JH3k4me3/target/$(H3k4me3_target)1.bed
	cp $(H3k4me3_target)2.bed $JH3k4me3/target/$(H3k4me3_target)2.bed
	bash $JJAMM.sh \
	    -g $(CHRSIZE) \
	    -o $JH3k4me3 \
	    -s $JH3k4me3/target \
	    -c $JH3k4me3/input

$JH3k36me3.bed: $(addsuffix .bed,$(H3k36me3_all)) $JJAMM.sh $(CHRSIZE)
	mkdir -p $JH3k36me3/input $JH3k36me3/target
	cp $(Uw_input).bed $JH3k36me3/input/$(Uw_input).bed
	cp $(H3k36me3_target)1.bed $JH3k36me3/target/$(H3k36me3_target)1.bed
	cp $(H3k36me3_target)2.bed $JH3k36me3/target/$(H3k36me3_target)2.bed
	bash $JJAMM.sh \
	    -g $(CHRSIZE) \
	    -o $JH3k36me3 \
	    -s $JH3k36me3/target \
	    -c $JH3k36me3/input \
	    -r region

$JPol2.bed: $(addsuffix .bed,$(Pol2_all)) $JJAMM.sh $(CHRSIZE)
	mkdir -p $JPol2/input $JPol2/target
	cp $(Pol2_input).bed $JPol2/input/$(Pol2_input).bed
	cp $(Pol2_target)1.bed $JPol2/target/$(Pol2_target)1.bed
	cp $(Pol2_target)2.bed $JPol2/target/$(Pol2_target)2.bed
	bash $JJAMM.sh \
	    -g $(CHRSIZE) \
	    -o $JPol2 \
	    -s $JPol2/target \
	    -c $JPol2/input

$MCtcf%.bed: $(Ctcf_target)%.sam $(Uw_input).sam
	macs2 callpeak \
	    -t $< \
	    -c $(Uw_input).sam \
	    -n $(subst .bed,,$@) \
	    --outdir $M \
	    -f SAM \
	    -g hs
	find $M -type f \
	    -not -name '*narrowPeak' \
	    -not -name '*broadPeak' | \
	xargs rm
	mv $(subst _peaks.narrowPeak,.bed,$@) $@


$MH3k4me3%.bed: $(H3k4me3_target)%.sam $(Uw_input).sam
	macs2 callpeak \
	    -t $< \
	    -c $(Uw_input).sam \
	    -n $(subst .bed,,$@) \
	    --outdir $M \
	    -f SAM \
	    -g hs
	find $M -type f \
	    -not -name '*narrowPeak' \
	    -not -name '*broadPeak' | \
	xargs rm
	mv $(subst _peaks.narrowPeak,.bed,$@) $@


$MH3k36me3%.bed: $(H3k36me3_target)%.sam $(Uw_input).sam
	macs2 callpeak \
	    -t $< \
	    -c $(Uw_input).sam \
	    -n $(subst .bed,,$@) \
	    --outdir $M \
	    -f SAM \
	    -g hs \
	    --broad
	find $M -type f \
	    -not -name '*narrowPeak' \
	    -not -name '*broadPeak' | \
	xargs rm
	mv $(subst _peaks.broadPeak,.bed,$@) $@


$MPol2%.bed: $(Pol2_target)%.sam $(Pol2_input).sam
	macs2 callpeak \
	    -t $< \
	    -c $(Pol2_input).sam \
	    -n $(subst .bed,,$@) \
	    --outdir $M \
	    -f SAM \
	    -g hs
	find $M -type f \
	    -not -name '*narrowPeak' \
	    -not -name '*broadPeak' | \
	xargs rm
	mv $(subst _peaks.narrowPeak,.bed,$@) $@


$ZCtcf.bed: $(addsuffix .map,$(Ctcf_all))
	$Zzerone -l \
	    --mock $(Uw_input).map \
	    --chip $(Ctcf_target)1.map \
	    --chip $(Ctcf_target)2.map \
	    > $@

$ZH3k4me3.bed: $(addsuffix .map,$(H3k4me3_all))
	$Zzerone -l \
	    --mock $(Uw_input).map \
	    --chip $(H3k4me3_target)1.map \
	    --chip $(H3k4me3_target)2.map \
	    > $@

$ZH3k36me3.bed: $(addsuffix .map,$(H3k36me3_all))
	$Zzerone -l \
	    --mock $(Uw_input).map \
	    --chip $(H3k36me3_target)1.map \
	    --chip $(H3k36me3_target)2.map \
	    > $@

$ZPol2.bed: $(addsuffix .map,$(Pol2_all))
	$Zzerone -l \
	    --mock $(Pol2_input).map \
	    --chip $(Pol2_target)1.map \
	    --chip $(Pol2_target)2.map \
	    > $@
