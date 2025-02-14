
# File: makefile
# Date: 14 February 2025
# Version 1: 
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Build EGGS.

CC?=gcc
CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/eggs: lib src/Main.o
	mkdir -p bin
	$(CC) $(LFLAGS) bin/eggs src/main.o lib/*.o -lz -lm 

src/Main.o:
	$(CC) $(CFLAGS) src/main.c -o src/main.o

.PHONY: lib
lib: 
	$(CC) $(CFLAGS) lib/bgzf.c -o lib/bgzf.o

.PHONY: clean
clean:
	rm lib/*.o src/Main.o bin/eggs