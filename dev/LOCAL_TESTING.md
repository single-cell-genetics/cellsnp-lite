# Local builds for development

Build out of tree in `build/` (already in `.gitignore`), so the source tree stays clean and
several builds (release, sanitiser, master) can sit side by side.

## macOS (Homebrew, Apple Silicon)

Verified with Homebrew htslib 1.23.1 and Apple clang.

```bash
brew install htslib xz

mkdir -p build/release && cd build/release
../../configure --with-htslib=/opt/homebrew \
    CPPFLAGS="-I/opt/homebrew/include" \
    LDFLAGS="-L/opt/homebrew/lib" \
    CFLAGS="-O2 -Wno-error=implicit-function-declaration" \
    LIBS="-ldeflate"
make -j8
./cellsnp-lite --version
```

Why each flag is needed:

- `CPPFLAGS` / `LDFLAGS`: `configure` checks for `lzma.h` and fails without the Homebrew prefix.
- `-Wno-error=implicit-function-declaration`: `src/thpool.c` calls `pthread_setname_np` without
  a visible declaration; recent Apple clang makes that an error. The call itself is correct on macOS.
- `LIBS="-ldeflate"`: `configure` links the static `libhts.a`, which Homebrew builds against
  libdeflate, so libdeflate has to be linked explicitly.

On Intel Macs replace `/opt/homebrew` with `/usr/local`.

## Linux

Not yet verified on this branch. The usual route is htslib built from source:

```bash
# htslib
git clone --recurse-submodules https://github.com/samtools/htslib.git
cd htslib && autoreconf -i && ./configure --prefix=$HOME/opt/htslib && make -j && make install

# cellsnp-lite
mkdir -p build/release && cd build/release
../../configure --with-htslib=$HOME/opt/htslib
make -j
```

`--with-htslib` needs an absolute path. If linking fails with undefined `libdeflate_*` symbols,
htslib was built with libdeflate: rerun `configure` with `LIBS="-ldeflate"`, as on macOS. If the
binary cannot find `libhts.so` at runtime, add `$HOME/opt/htslib/lib` to `LD_LIBRARY_PATH`.

## Sanitiser builds

Same `configure` line, different `make` flags. Verified on macOS.

```bash
# AddressSanitizer
make -j8 CFLAGS="-O1 -g -fsanitize=address -fno-omit-frame-pointer -Wno-error=implicit-function-declaration" \
    LDFLAGS="-L/opt/homebrew/lib -fsanitize=address"

# ThreadSanitizer
make -j8 CFLAGS="-O1 -g -fsanitize=thread -Wno-error=implicit-function-declaration" \
    LDFLAGS="-L/opt/homebrew/lib -fsanitize=thread"
```

Run `make clean` first when switching an existing build directory between flavours.

TSan on master already reports three races in the vendored `src/thpool.c` (pool start-up and
shutdown) and aborts at exit. These are pre-existing; anything outside `thpool.c` is new.

## Comparing against master

Performance changes should not change output. Build master next to the branch and byte-compare:

```bash
git worktree add ../cellsnp-lite-master master
mkdir -p build/master && cd build/master
../../../cellsnp-lite-master/configure <same flags as above>
make -j8
```

Run both binaries with the same arguments into different `-O` directories and `cmp` every output
file (`zcat` first with `--gzip`). Cover at least `-R`, `-T`, mode 2, `--genotype`, `-p 1` and
`-p > 1`, `--UMItag None` and multiple samples with `-I`.

The 10x test data used upstream (130 MB) has a BAM, 400 barcodes and a hg19 SNP subset:
<https://sourceforge.net/projects/cellsnp/files/testdata/cellSNP_testdata_10x.zip/download>
