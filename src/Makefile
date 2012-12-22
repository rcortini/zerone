OBJECTS = fwdb.o viterbi.o HummingBee.so 

all: HummingBee.so

clean:
	rm -f $(OBJECTS)

HummingBee.so: $(OBJECTS)
	R CMD SHLIB fwdb.o viterbi.o -o HummingBee.so -lm
	#gcc -std=gnu99 -shared -L/usr/local/lib64 -o HummingBee.so \
	#	ffwdb.o fwdb.o viterbi.o viterbi.o -lgfortran -lm \
	#	-L/soft/general/R-2.15.0/lib64/R/lib -lR

fwdb.o: fwdb.c
	#gcc -std=gnu99 -fPIC -g -O3 -c fwdb.c -lpthread -Wall
	gcc -std=gnu99 -fPIC -g -O0 -c fwdb.c -lpthread -Wall

viterbi.o: viterbi.c
	gcc -std=gnu99 -fPIC -g -O3 -c viterbi.c -lpthread -Wall
