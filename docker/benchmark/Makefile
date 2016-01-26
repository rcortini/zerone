BW2BURL= http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/
GRO= GSM1480325_K562_GROseq_plus.bed GSM1480325_K562_GROseq_minus.bed

all: $(GRO)
	#make -C ../discretize

bigWigToBedGraph:
	wget $(BW2BURL)$@
	chmod +x $@

GSM1480325_K562_GROseq_%.bed: GSM1480325_K562_GROseq_%.bigWig bigWigToBedGraph
	./bigWigToBedGraph $< $@
