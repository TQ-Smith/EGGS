
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
	$(CC) $(CFLAGS) src/Missingness.c -o src/Missingness.o

.PHONY: lib 
lib: lib/kstring.o lib/gsl lib/fftpack/fft.o

lib/kstring.o: lib/malloc_wrap.o
	$(CC) $(CFLAGS) lib/kstring.c -o lib/kstring.o

lib/malloc_wrap.o:
	$(CC) $(CFLAGS) lib/malloc_wrap.c -o lib/malloc_wrap.o

.PHONY: lib/gsl
lib/gsl:
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/beta.c -o lib/gsl/beta.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/gamma.c -o lib/gsl/gamma.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/gausszig.c -o lib/gsl/gausszig.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/error.c -o lib/gsl/error.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/message.c -o lib/gsl/message.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/stream.c -o lib/gsl/stream.o

lib/fftpack/fft.o:
	$(CC) $(CFLAGS) lib/fftpack/fft.c -o lib/fftpack/fft.o

.PHONY: clean
clean:
	rm lib/*.o src/*.o lib/gsl/*.o lib/fftpack/*.o bin/eggs