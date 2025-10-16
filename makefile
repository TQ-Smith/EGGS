
# File: makefile
# Date: 13 May 2025
# Version 1: 
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Build EGGS.

CC?=gcc
CFLAGS = -c -Wall -g -I lib
LFLAGS = -g -o

bin/eggs: lib src
	mkdir -p bin
	$(CC) $(LFLAGS) bin/eggs src/*.o lib/*.o lib/gsl/*.o -lz -lm

.PHONY: src
src: src/Main.o

src/Main.o: src/Interace.o src/GenotypeParser.o src/Missingness.o
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

src/Interace.o:
	$(CC) $(CFLAGS) src/Interface.c -o src/Interace.o

src/GenotypeParser.o:
	$(CC) $(CFLAGS) src/GenotypeParser.c -o src/GenotypeParser.o

src/Missingness.o:
	$(CC) $(CFLAGS) -DHAVE_INLINE src/Missingness.c -o src/Missingness.o

.PHONY: lib 
lib: lib/kstring.o lib/gsl

lib/kstring.o: 
	$(CC) $(CFLAGS) lib/kstring.c -o lib/kstring.o

.PHONY: lib/gsl
lib/gsl:
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/error.c -o lib/gsl/error.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/message.c -o lib/gsl/message.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/stream.c -o lib/gsl/stream.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/default.c -o lib/gsl/default.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/rng.c -o lib/gsl/rng.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/mt.c -o lib/gsl/mt.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/types.c -o lib/gsl/types.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/gauss.c -o lib/gsl/gauss.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/beta.c -o lib/gsl/beta.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/gamma.c -o lib/gsl/gamma.o

.PHONY: clean
clean:
	rm lib/*.o src/*.o lib/gsl/*.o
