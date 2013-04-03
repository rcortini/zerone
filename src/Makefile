OBJECTS= estep.o viterbi.o indexts.o sha1.o loglik.o HummingBee.so 

all: HummingBee.so

clean:
	rm -f $(OBJECTS)

HummingBee.so: $(OBJECTS)
	R CMD SHLIB estep.o viterbi.o indexts.o	loglik.o \
	sha1.o -o HummingBee.so -lpthread -lm

loglik.o: loglik.c loglik.h
	gcc -std=gnu99 -fPIC -g -O3 -c loglik.c -Wall

estep.o: estep.c estep.h
	gcc -std=gnu99 -fPIC -g -O3 -c estep.c -Wall

viterbi.o: viterbi.c
	gcc -std=gnu99 -fPIC -g -O3 -c viterbi.c -Wall

indexts.o: indexts.c indexts.h
	gcc -std=gnu99 -fPIC -g -O3 -c indexts.c -Wall

sha1.o: sha1.c sha1.h
	gcc -std=gnu99 -fPIC -g -O3 -c sha1.c -Wall
