vpath %.c ..
vpath %.c lib
vpath %.h ..
vpath %.h lib

P= runtests

OBJECTS= libunittest.so xxhash.o sam.o bgzf.o hfile.o \
	 tests_zinm.o tests_hmm.o tests_predict.o tests_zerone.o \
	 tests_parse.o tests_utils.o
SOURCES= zinm.c hmm.c predict.c zerone.c parse.c utils.c

CC= gcc
INCLUDES= -I.. -Ilib
COVERAGE= -fprofile-arcs -ftest-coverage
CFLAGS= -std=gnu99 -g -Wall -O0 $(INCLUDES) $(COVERAGE)
LDLIBS= -L`pwd` -Wl,-rpath=`pwd` -lunittest -lz -lm

$(P): $(OBJECTS) $(SOURCES) runtests.c
	$(CC) $(CFLAGS) runtests.c $(OBJECTS) $(LDLIBS) -o $@

clean:
	rm -f $(P) $(OBJECTS) *.gcda *.gcno *.gcov gmon.out .inspect.gdb

libunittest.so: unittest.c
	$(CC) -fPIC -shared $(CFLAGS) -o libunittest.so lib/unittest.c

test: $(P)
	./$(P)

inspect: $(P)
	gdb --command=.inspect.gdb --args $(P)

debug: $(P)
	gdb --command=.debug.gdb --args $(P)

valgrind: $(P)
	#valgrind --vgdb=yes --vgdb-error=0 ./$(P) --debug
	valgrind --leak-check=full --show-leak-kinds=all ./$(P) --debug
