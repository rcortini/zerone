vpath %.c ..
vpath %.h ..

P=test
OBJECTS= pso.o mnmultinom.o hmm.o utils.o sha1.o
CFLAGS= -I.. `pkg-config --cflags glib-2.0` -g -pg -Wall -std=gnu99 \
	          -O3 -fprofile-arcs -ftest-coverage
LDLIBS= `pkg-config --libs glib-2.0` -lpthread -lm
CC=gcc
$(P): $(OBJECTS)

clean:
	rm -f test *.o *.gcda *.gcno *.gcov gmon.out analysis.txt \
		callgrind.out.* cache.txt

prof:
	./test -p && gprof -F mnmultinom_prob -F block_fwdb test \
		gmon.out > analysis.txt
	gcov *.gcda
