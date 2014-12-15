SRC_DIR= ../../src
INC_DIR= ../../src

SOURCES= $(addprefix $(SRC_DIR)/,$(SOURCE_FILES))
INCLUDES= $(addprefix -I, $(INC_DIR))
OBJECTS= Rjahmm.o jahmm.o zinb.o hmm.o utils.o xxhash.o

CC= gcc
CFLAGS= -std=gnu99 -fPIC -g -Wall -O3

all: Rjahmm.so

%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

Rjahmm.o: Rjahmm.c
	R CMD SHLIB $(SRC_DIR)/Rjahmm.c -o Rjahmm.o

Rjahmm.so: $(OBJECTS) $(SOURCES)
	R CMD SHLIB $(SOURCES) $(OBJECTS) -o Rjahmm.so

clean:
	rm -f $(OBJECTS) Rjahmm.so
