OBJECTS = estep.o viterbi.o loglik.o HummingBee.so 

all: HummingBee.so

clean:
	rm -f $(OBJECTS)

HummingBee.so: $(OBJECTS)
	R CMD SHLIB estep.o viterbi.o loglik.o -o HummingBee.so -lm

loglik.o: loglik.c
	gcc -std=gnu99 -fPIC -g -O3 -c estep.c -Wall -lm

estep.o: estep.c
	gcc -std=gnu99 -fPIC -g -O3 -c estep.c -Wall -lm

viterbi.o: viterbi.c
	gcc -std=gnu99 -fPIC -g -O3 -c viterbi.c -Wall -lm
