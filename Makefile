CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
BINDIR=/usr/local/bin

htslib_dir=../htslib
htslib_include_dir=$(htslib_dir)
htslib_lib_dir=$(htslib_dir)
scripts=cellsnp.c thpool.c
headers=general_util.h cellsnp_util.h thpool.h kvec.h

all: cellSNP

cellSNP: $(scripts) $(headers)
	$(CC) $(CFLAGS) $(scripts) -o $@ -I $(htslib_include_dir) -L $(htslib_lib_dir) -lz -lm -lhts -pthread

install: all
	install cellSNP $(BINDIR)

clean:
	-rm *.o a.out cellSNP
