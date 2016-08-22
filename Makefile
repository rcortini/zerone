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

P= zerone

SRC_DIR= src
INC_DIR= src

OBJECT_FILES= bgzf.o sam.o hfile.o hmm.o utils.o xxhash.o zerone.o \
      zinm.o parse.o snippets.o
SOURCE_FILES= main.c predict.c

OBJECTS= $(addprefix $(SRC_DIR)/,$(OBJECT_FILES))
SOURCES= $(addprefix $(SRC_DIR)/,$(SOURCE_FILES))
INCLUDES= $(addprefix -I, $(INC_DIR))

CC= gcc
CFLAGS= -std=gnu99 -Wall
LDLIBS= -lm

all: CFLAGS += -O3
all: $(P)

debug: CFLAGS += -DDEBUG -g -O0
debug: $(P)

atac: CFLAGS += -O3 -DATAC_SEQ
atac: $(P)_atac

$(P): $(OBJECTS) $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@ -lz -lm

$(P)_atac: $(OBJECTS) $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@ -lz -lm


$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(P)

