CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
BINDIR=/usr/local/bin

scripts=cellsnp.c thpool.c
headers=general_util.h cellsnp_util.h thpool.h kvec.h

all: cellSNP

cellSNP: $(scripts) $(headers)
	$(CC) $(CFLAGS) $(scripts) -o $@ -lz -lm -lhts -pthread

install: all
	install cellSNP $(BINDIR)

clean:
	-rm *.o a.out cellSNP