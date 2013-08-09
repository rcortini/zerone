vpath %.c ..
vpath %.h ..

P= testset
OBJECTS= pso.o mnmultinom.o hmm.o utils.o xxhash.o
CFLAGS= -I.. `pkg-config --cflags glib-2.0` -g -pg -Wall -std=gnu99 \
	          -O3 -fstrict-aliasing -fprofile-arcs -ftest-coverage
LDLIBS= `pkg-config --libs glib-2.0` -lpthread -lm
CC= gcc
$(P): $(OBJECTS)

clean:
	rm -f testset *.o *.gcda *.gcno *.gcov gmon.out analysis.txt \
		callgrind.out.* cache.txt

test:
	gtester --verbose --keep-going testset

testperf:
	gtester --verbose --keep-going -m=perf testset
	gprof -F mnmultinom_prob -F block_fwdb testset \
		gmon.out > analysis.txt
	gcov *.gcda
