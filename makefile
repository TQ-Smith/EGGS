
# File: makefile
# Date: 14 February 2025
# Version 1: 
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Build EGGS.

CC?=gcc
CFLAGS = -c -Wall -g -DUSE_MALLOC_WRAPPERS -I lib
LFLAGS = -g -o

bin/eggs: lib src
	mkdir -p bin
	$(CC) $(LFLAGS) bin/eggs src/*.o lib/*.o -lz -lm 

.PHONY: src
src: src/Main.o

src/Main.o: src/Interace.o src/GenotypeParser.o 
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

src/Interace.o:
	$(CC) $(CFLAGS) src/Interface.c -o src/Interace.o

src/GenotypeParser.o:
	$(CC) $(CFLAGS) src/GenotypeParser.c -o src/GenotypeParser.o

.PHONY: lib 
lib: lib/kstring.o

lib/kstring.o: lib/malloc_wrap.o
	$(CC) $(CFLAGS) lib/kstring.c -o lib/kstring.o

lib/malloc_wrap.o:
	$(CC) $(CFLAGS) lib/malloc_wrap.c -o lib/malloc_wrap.o

.PHONY: clean
clean:
	rm lib/*.o src/*.o bin/eggs