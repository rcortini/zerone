SRC_DIR= ../../src
INC_DIR= ../../src

INCLUDES= $(addprefix -I, $(INC_DIR))
OBJECTS= zerone.o zinb.o hmm.o utils.o xxhash.o

CC= gcc
CFLAGS= -std=gnu99 -fPIC -g -Wall -O3

all: Rzerone.so

%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

Rzerone.so: $(OBJECTS)
	R CMD SHLIB $(INCLUDES) Rzerone.c $(OBJECTS)

clean:
	rm -f $(OBJECTS) Rzerone.o Rzerone.so
