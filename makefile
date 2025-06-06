
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
	$(CC) $(LFLAGS) bin/eggs src/*.o lib/*.o lib/gsl/*.o lib/kissfft/*.o -lz -lm -lpthread

.PHONY: src
src: src/Main.o

src/Main.o: src/Interace.o src/GenotypeParser.o src/Missingness.o
	$(CC) $(CFLAGS) -DUSE_MALLOC_WRAPPERS src/Main.c -o src/Main.o

src/Interace.o:
	$(CC) $(CFLAGS) -DUSE_MALLOC_WRAPPERS src/Interface.c -o src/Interace.o

src/GenotypeParser.o:
	$(CC) $(CFLAGS) -DUSE_MALLOC_WRAPPERS src/GenotypeParser.c -o src/GenotypeParser.o

src/Missingness.o: src/Sort.o
	$(CC) $(CFLAGS) -DUSE_MALLOC_WRAPPERS -DHAVE_INLINE src/Missingness.c -o src/Missingness.o

src/Sort.o:
	$(CC) $(CFLAGS) src/Sort.c -o src/Sort.o

.PHONY: lib 
lib: lib/kstring.o lib/gsl lib/kissfft

lib/kstring.o: lib/malloc_wrap.o
	$(CC) $(CFLAGS) -DUSE_MALLOC_WRAPPERS lib/kstring.c -o lib/kstring.o

lib/malloc_wrap.o:
	$(CC) $(CFLAGS) lib/malloc_wrap.c -o lib/malloc_wrap.o

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
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/gausszig.c -o lib/gsl/gausszig.o

.PHONY: lib/kissfft
lib/kissfft:
	$(CC) $(CFLAGS) lib/kissfft/kfc.c -o lib/kissfft/kfc.o
	$(CC) $(CFLAGS) lib/kissfft/kiss_fft.c -o lib/kissfft/kiss_fft.o
	$(CC) $(CFLAGS) lib/kissfft/kiss_fftnd.c -o lib/kissfft/kiss_fftnd.o
	$(CC) $(CFLAGS) lib/kissfft/kiss_fftndr.c -o lib/kissfft/kiss_fftndr.o
	$(CC) $(CFLAGS) lib/kissfft/kiss_fftr.c -o lib/kissfft/kiss_fftr.o

.PHONY: clean
clean:
	rm lib/*.o src/*.o lib/gsl/*.o lib/kissfft/*.o bin/eggs