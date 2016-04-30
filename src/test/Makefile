# Copyright 2015, 2016 Pol Cusco and Guillaume Filion
#
# This file is part of Zerone.
#
# Zerone is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Zerone is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Zerone. If not, see <http://www.gnu.org/licenses/>.

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
LDLIBS= -L`pwd` -Wl,-rpath=`pwd` -lunittest -lz -lm -lpthread

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
