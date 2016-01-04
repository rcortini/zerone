P= zerone

SRC_DIR= src
INC_DIR= src

OBJECT_FILES= bgzf.o sam.o hfile.o hmm.o utils.o xxhash.o zerone.o \
      zinm.o parse.o
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

$(P): $(OBJECTS) $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@ -lz -lm

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(P)

