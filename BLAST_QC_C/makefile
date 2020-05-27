CC=gcc
CFLAGS=-I./libxml2
LFLAGS=-lxml2
RM=rm -rf
OUT=BLASTQC

all: build

build: BLASTQC.o vector.o
	$(CC) $(CFLAGS) -o $(OUT) BLASTQC.c vector.c $(LFLAGS)
	$(RM) *.o

debug: 
	CFLAGS+=-DDEBUG_ON

debug: build

BLASTQC.o: BLASTQC.c vector.h BLASTQC.h
	$(CC) $(CFLAGS) -c BLASTQC.c $(LFLAGS)

vector.o: vector.c vector.h
	$(CC) $(CFLAGS) -c vector.c $(LFLAGS)

clean:
	$(RM) *.o $(OUT)
