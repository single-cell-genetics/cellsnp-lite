/* File operatoins API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef SZ_JFILE_H
#define SZ_JFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include "htslib/kstring.h"        // do not use "kstring.h" as it's different from "htslib/kstring.h"
#include "htslib/bgzf.h"
#include "config.h"

/*
 * File structure with simple output buffer (Support bgzip)
 */

/* let gzgets() to be compatible with kgetline(). */
static inline char* jf_gzgets(char *buf, int len, gzFile fp) { return gzgets(fp, buf, len); }
static inline int jf_bgzf_gets(BGZF *fp, kstring_t *s) {
    return bgzf_getline(fp, '\n', s) <= -1 ? -1 : 0;
}

#define JF_GZIP 1
#define JF_BGZIP 2

#ifndef JF_ZIP_TYPE
#define JF_ZIP_TYPE JF_BGZIP
#endif

#if (JF_ZIP_TYPE == JF_BGZIP)
    #define JF_ZFILE BGZF*
    #define jf_zopen(fn, fm) bgzf_open(fn, fm)
    #define jf_zclose(zfp) bgzf_close(zfp)
    #define jf_zread(zfp, buf, len) bgzf_read(zfp, buf, len)
    #define jf_zgetln(zfp, ks) jf_bgzf_gets(zfp, ks)
    #define jf_zwrite(zfp, buf, len) bgzf_write(zfp, buf, len)
#else
    #define JF_ZFILE gzFile
    #define jf_zopen(fn, fm) gzopen(fn, fm)
    #define jf_zclose(zfp) gzclose(zfp)
    #define jf_zread(zfp, buf, len) gzread(zfp, buf, len)
    #define jf_zgetln(zfp, ks) kgetline(ks, (kgets_func*) jf_gzgets, zfp)
    #define jf_zwrite(zfp, buf, len) gzwrite(zfp, buf, len)
#endif

/*@abstract    Structure used for outputed files with Simple Output Buffer.
@param fn      Filename.
@param fm      Filemode.
@param zfp     File pointer for zipped mode.
@param fp      File pointer for non-zipped mode.
@param is_zip  If the outputed file should be zipped.
@param is_tmp  If the outputed file is a tmp file. It's only used by mtx file and has no effect on I/O.
@param is_open If the outputed file is open.
@param buf     Mimic Output Buffer.
@param bufsize Size of buffer.
@note          1. The @p fn should be valid pointer coming from strdup().
               2. The @p fm points to const string, so do not free it!
               3. Output buffer is inside the structure.
@TODO  Add is_error to save the state that if I/O has error.
 */
typedef struct {
    char *fn;
    char *fm;
    JF_ZFILE zfp;
    FILE *fp;
    uint8_t is_zip, is_tmp, is_open;
    kstring_t ks, *buf;
    size_t bufsize;
} jfile_t;

/*@abstract  Initialize the jfile_t structure.
@param       Void.
@return      Pointer of jfile_t if success, NULL otherwise.
 */
static inline jfile_t* jf_init(void) { 
    jfile_t *p = (jfile_t*) calloc(1, sizeof(jfile_t));
    if (p) { p->buf = &p->ks; p->bufsize = 1048576; }   // p->ks has been initialized after calling calloc(). Double initializing will cause error.
    return p;
}

static inline void jf_destroy(jfile_t* p) {
    if (p) {
        if (p->is_open) {
            if (p->is_zip) { jf_zclose(p->zfp); p->zfp = NULL; }
            else { fclose(p->fp); p->fp = NULL; }
        }
        ks_free(p->buf);
        free(p->fn); free(p);
    }
}

static inline void jf_set_bufsize(jfile_t *p, size_t bufsize) { p->bufsize = bufsize; }

/* If the jfile_t is open. 0:no; 1:yes. */
static inline int jf_isopen(jfile_t *p) { return p->is_open; }

/*@abstract  Open jfile_t for reading or writing.
@param p     Pointer of jfile_t.
@param mode  Open mode as in fopen(). if NULL, use default mode inside jfile_t.
@return      1 if success, 0 if already open, -1 if error. 
 */
static inline int jf_open(jfile_t *p, char *mode) {
    if (p->is_open) { return 0; }
    char *fm = mode ? mode : p->fm;
    if (p->is_zip) {
        if (NULL == (p->zfp = jf_zopen(p->fn, fm))) { return -1; }
        else { p->is_open = 1; return 1; }
    } else if (NULL == (p->fp = fopen(p->fn, fm))) {
        return -1;
    } else { p->is_open = 1; return 1; }
}

/*@abstract  Read from jfile_t without Input buffer.
@param p     Pointer of jfile_t.
@param buf   Buffer where the read content will be pushed into.
@param len   Size of content to be read.
@return      Value returned as gzread/fread. */
static inline ssize_t jf_read(jfile_t *p, char *buf, size_t len) {
    return p->is_zip ? jf_zread(p->zfp, buf, len) : fread(buf, 1, len, p->fp);
}

/*@abstract  Get a line from jfile_t without Input buffer.
@param p     Pointer of jfile_t.
@param s     Pointer of kstring_t which saves the line.
@return      Value returned as kgetline(). 0 if success, EOF if end-of-file or error.

@note        kgetline() is defined in htslib/kstring.h 
*/
static inline int jf_getln(jfile_t *p, kstring_t *s) {
    return p->is_zip ? jf_zgetln(p->zfp, s) : kgetline(s, (kgets_func*) fgets, p->fp);
}

/*@abstract  Output functions below are like kxxx() functions in htslib/kstring.h with Output buffer.
@return      Value returned are like in fxxx()/kxxx() functions (fputc()/kputc() etc.), which are 
             the real size outputed if success, EOF if error. */

/*@abstract  Flush buffer to stream.
@param p     Pointer of jfile_t.
@return      0 if success, EOF otherwise.
 */
static inline int jf_flush(jfile_t *p) {
    ssize_t l, l0 = ks_len(p->buf);
    l = p->is_zip ? jf_zwrite(p->zfp, ks_str(p->buf), ks_len(p->buf)) : fwrite(ks_str(p->buf), 1, ks_len(p->buf), p->fp);
    ks_clear(p->buf);
    return l == l0 ? 0 : EOF;
}

static inline int jf_printf(jfile_t *p, const char *fmt, ...) {
    va_list ap;
    int l;
    va_start(ap, fmt);
    l = kvsprintf(p->buf, fmt, ap);
    va_end(ap);
    if (ks_len(p->buf) >= p->bufsize && jf_flush(p) < 0) { return EOF; }
    return l;
}

static inline int jf_putc(int c, jfile_t *p) { 
    int l;
    l = kputc(c, p->buf);
    if (ks_len(p->buf) >= p->bufsize && jf_flush(p) < 0) { return EOF; }
    return l;
}

static inline int jf_putc_(int c, jfile_t *p) { 
    int l;
    l = kputc_(c, p->buf);
    if (ks_len(p->buf) >= p->bufsize && jf_flush(p) < 0) { return EOF; }
    return l;
}

static inline int jf_puts(const char *s, jfile_t *p) { 
    int l;
    l = kputs(s, p->buf); 
    if (ks_len(p->buf) >= p->bufsize && jf_flush(p) < 0) { return EOF; }
    return l;
}

static inline int jf_write(jfile_t *p, char *buf, size_t len) { 
    int l; 
    l = kputsn(buf, len, p->buf);
    if (ks_len(p->buf) >= p->bufsize && jf_flush(p) < 0) { return EOF; }
    return l;
}

/*@abstract  Close jfile_t but does not destroy it.
@param p     Pointer of jfile_t.
@return      0 if success, EOF otherwise.

@note        Even fail, the jfile_t will still be set to not open.
 */
static inline int jf_close(jfile_t *p) {
    int ret = 0;
    if (p->is_open) {
        if (ks_len(p->buf) && jf_flush(p) < 0) { ret = EOF; } // only for write mode.
        if (p->is_zip) { jf_zclose(p->zfp); p->zfp = NULL; }
        else { fclose(p->fp); p->fp = NULL; }
        p->is_open = 0;
    } 
    return ret;
}

/*@abstract   Remove file.
@param p      Pointer of jfile_t to be removed.
@return       1 if success, 0 if file does not exist, -1 if failure.
 */
static inline int jf_remove(jfile_t *p) {
    if (0 != access(p->fn, F_OK)) { return 0; }
    if (remove(p->fn) < 0) { return -1; }
    return 1;
}

/*@abstract   Remove all files in array.
@param fs     Pointer of array of jfile_t to be removed.
@param n      Size of array.
@return       Num of files have been removed if no error, -1 otherwise.
 */
static inline int jf_remove_all(jfile_t **fs, const int n) {
    int i, j, ret;
    for (i = 0, j = 0; i < n; i++) {
        if ((ret = jf_remove(fs[i])) < 0) { return -1; }
        else { j += ret; }
    }
    return j;
}

/* 
* File Functions
 */

/*@abstract   Join together two pathes.
@param p1     Pointer to the first path.
@param p2     Pointer to the second path.
@return       Pointer to the joined path if success, NULL otherwise.

@note         Only works for Unix system as the path seperator used in this function is '/'.
 */
static inline char* join_path(const char *p1, const char *p2) {
    kstring_t ks = KS_INITIALIZE, *s = &ks;
    char *p = NULL;
    int n1;
    if (NULL == p1 || (n1 = strlen(p1)) <= 0 || NULL == p2) { ks_free(s); return NULL; }
    kputs(p1, s);
    if (p1[n1 - 1] != '/') { kputc('/', s); }
    kputs(p2, s); 
    p = strdup(ks_str(s));
    ks_free(s);
    return p;
}

#define TMP_BUFSIZE 1048576
/*@abstract      Merge several files into one.
@param in_fn     Names of input files to be merged.
@param n         Num of input files.
@param out_fn    Name of output file.
@return          Num of files have been merged.
 */
static int merge_files(char **in_fn, const int n, const char *out_fn) {
    char buf[TMP_BUFSIZE];
    FILE *in = NULL, *out = NULL;
    int i = 0;
    size_t m;
    if (NULL == in_fn || NULL == out_fn) { return 0; }
    if (NULL == (out = fopen(out_fn, "w"))) { fprintf(stderr, "[E::%s] could not open '%s'\n", __func__, out_fn); return 0; }
    for (i = 0; i < n; i++) {
        if (NULL == (in = fopen(in_fn[i], "r"))) { fprintf(stderr, "[E::%s] could not open '%s'\n", __func__, in_fn[i]); goto fail; }
        while ((m = fread(buf, 1, TMP_BUFSIZE, in)) > 0) { fwrite(buf, 1, m, out); }
        if (ferror(in)) { goto fail; }
        else { fclose(in); in = NULL; }
    }
    fclose(out);
    return i;
  fail:
    if (in) fclose(in);
    if (out) fclose(out);
    return i;
}
#undef TMP_BUFSIZE

/*@abstract  Remove file.
@param fn    Name of file to be removed.
@return      1 if success, 0 if file does not exist, -1 if failure.
 */
static inline int remove_file(char *fn) {
    if (0 != access(fn, F_OK)) { return 0; }
    if (remove(fn) < 0) { return -1; }
    return 1;
}

/*@abstract   Remove several files.
@param fn     Names of files to be removed.
@param n      Num of files.
@return       Num of files have been removed if no error, -1 otherwise.
 */
static inline int remove_files(char **fn, const int n) {
    int i, j, ret;
    for (i = 0, j = 0; i < n; i++) {
        if ((ret = remove_file(fn[i])) < 0) { return -1; }
        else { j += ret; }
    }
    return j;
}

#endif
