
# File: makefile
# Date: 14 February 2025
# Version 1: 
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Build EGGS.

CC?=gcc
CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/eggs: src/main.o
	mkdir -p bin
	$(CC) $(LFLAGS) bin/eggs src/main.o -lz -lm 

src/main.o:
	$(CC) $(CFLAGS) src/main.c -o src/main.o

.PHONY: clean
clean:
	rm src/main.o bin/eggs