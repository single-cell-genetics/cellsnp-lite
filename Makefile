htslib_dir=../htslib
htslib_include_dir=$(htslib_dir)
htslib_lib_dir=$(htslib_dir)

CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function -fgnu89-inline -I$(htslib_include_dir)
LDFLAGS=-L$(htslib_lib_dir)

BIN_DIR=/usr/local/bin
BIN_NAME=cellsnp-lite

src_dir=src
scripts=$(src_dir)/cellsnp.c $(src_dir)/csp_fetch.c $(src_dir)/csp_pileup.c $(src_dir)/csp.c $(src_dir)/jfile.c $(src_dir)/jsam.c $(src_dir)/jstring.c $(src_dir)/mplp.c $(src_dir)/snp.c $(src_dir)/thpool.c
headers=$(src_dir)/config.h $(src_dir)/csp.h $(src_dir)/jfile.h $(src_dir)/jmempool.h $(src_dir)/jnumeric.h $(src_dir)/jsam.h $(src_dir)/jstring.h $(src_dir)/kvec.h $(src_dir)/mplp.h $(src_dir)/snp.h $(src_dir)/thpool.h

all: $(BIN_NAME)

$(BIN_NAME): $(scripts) $(headers)
	$(CC) $(CFLAGS) $(LDFLAGS) $(scripts) -o $@ -lz -lm -lhts -pthread

install: all
	install $(BIN_NAME) $(BIN_DIR)

clean:
	-rm -f *.o a.out $(BIN_NAME)
