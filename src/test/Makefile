vpath %.c ..
vpath %.h ..

P=test
OBJECTS= loglik.o indexts.o sha1.o
CFLAGS= -I.. `pkg-config --cflags glib-2.0` -g -pg -Wall -std=gnu99 \
	          -O3 -fprofile-arcs -ftest-coverage
LDLIBS= `pkg-config --libs glib-2.0` -lpthread -lm
CC=gcc
$(P): $(OBJECTS)

clean:
	rm -f test *.o *.gcda *.gcno *.gcov gmon.out
