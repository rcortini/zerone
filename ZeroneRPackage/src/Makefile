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

SRC_DIR= ../../src
INC_DIR= ../../src

INCLUDES= $(addprefix -I, $(INC_DIR))
OBJECTS= zerone.o zinm.o hmm.o utils.o xxhash.o

CC= gcc
CFLAGS= -std=gnu99 -fPIC -Wall

all: CFLAGS += -O3
all: Rzerone.so

debug: CFLAGS += -DDEBUG -g -O0
debug: SHLIBFLAGS += -d
debug: Rzerone.so

%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

Rzerone.so: $(OBJECTS)
	R CMD SHLIB $(SHLIBFLAGS) $(INCLUDES) Rzerone.c $(OBJECTS)

clean:
	rm -f $(OBJECTS) Rzerone.o Rzerone.so
