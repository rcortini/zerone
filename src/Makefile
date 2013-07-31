OBJECTS= pso.o mnmultinom.o hmm.o utils.o xxhash.o
CC= gcc
CFLAGS= -std=gnu99 -fPIC -g -O3 -Wall

all: jahmm.so

clean:
	rm -f $(OBJECTS) jahmm.so

jahmm.so: $(OBJECTS)
	R CMD SHLIB pso.o mnmultinom.o hmm.o utils.o xxhash.o \
	-o jahmm.so -lpthread -lm
