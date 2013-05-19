OBJECTS= pso.o mnmultinom.o hmm.o utils.o sha1.o
CC= gcc
CFLAGS= -std=gnu99 -fPIC -g -O3 -Wall

all: yahmm.so

clean:
	rm -f $(OBJECTS) yahmm.so

yahmm.so: $(OBJECTS)
	R CMD SHLIB pso.o mnmultinom.o hmm.o utils.o sha1.o \
	-o yahmm.so -lpthread -lm
