htslib_dir=../htslib
htslib_include_dir=$(htslib_dir)
htslib_lib_dir=$(htslib_dir)

CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function -I$(htslib_include_dir)
LDFLAGS=-L$(htslib_lib_dir)

BIN_DIR=/usr/local/bin
BIN_NAME=cellsnp-lite

scripts=cellsnp.c thpool.c
headers=general_util.h cellsnp_util.h thpool.h kvec.h

all: $(BIN_NAME)

$(BIN_NAME): $(scripts) $(headers)
	$(CC) $(CFLAGS) $(LDFLAGS) $(scripts) -o $@ -lz -lm -lhts -pthread

install: all
	install $(BIN_NAME) $(BIN_DIR)

clean:
	-rm -f *.o a.out $(BIN_NAME)
