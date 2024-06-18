// jfile.h - file operations.


#ifndef SZ_JFILE_H
#define SZ_JFILE_H

#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <zlib.h>
#include "htslib/kstring.h"        // do not use "kstring.h" as it's different from "htslib/kstring.h"
#include "htslib/bgzf.h"
#include "config.h"

/*
 * File structure with simple output buffer (Support bgzip)
 */

/* let gzgets() to be compatible with kgetline(). */
char* jf_gzgets(char *buf, int len, gzFile fp);
int jf_bgzf_gets(BGZF *fp, kstring_t *s);

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

// Structure used for outputed files with Simple Output Buffer.
// @TODO  Add is_error to save the state that if I/O has error.
typedef struct {
    char *fn;       // File name. should come from strdup().
    char *fm;       // File mode. point to const string, do not free it!
    JF_ZFILE zfp;   // File pointer for zipped mode.
    FILE *fp;       // File pointer for non-zipped mode.
    uint8_t is_zip, is_tmp, is_open;  // is_*: for output file. is_tmp: only for mtx file.
    kstring_t ks, *buf;
    size_t bufsize;
} jfile_t;

jfile_t* jf_init(void); 
void jf_destroy(jfile_t* p);
void jf_set_bufsize(jfile_t *p, size_t bufsize);
int jf_isopen(jfile_t *p);

/*!@func
@abstract    Open jfile_t for reading or writing.
@param p     Pointer of jfile_t.
@param mode  Open mode as in fopen(). if NULL, use default mode inside jfile_t.
@return      1 if success, 0 if already open, -1 if error. 
 */
int jf_open(jfile_t *p, char *mode);

/*!@func
@abstract    Read from jfile_t without Input buffer.
@param p     Pointer of jfile_t.
@param buf   Buffer where the read content will be pushed into.
@param len   Size of content to be read.
@return      Value returned as gzread/fread. 
*/
ssize_t jf_read(jfile_t *p, char *buf, size_t len);

/*!@func
@abstract  Get a line from jfile_t without Input buffer.
@param p   Pointer of jfile_t.
@param s   Pointer of kstring_t which saves the line.
@return    Value returned as kgetline(). 0 if success, EOF if end-of-file or error.
@note      kgetline() is defined in htslib/kstring.h 
*/
int jf_getln(jfile_t *p, kstring_t *s);

/*@abstract  Output functions below are like kxxx() functions in htslib/kstring.h with Output buffer.
@return      Value returned are like in fxxx()/kxxx() functions (fputc()/kputc() etc.), which are 
             the real size outputed if success, EOF if error. */

/*!@func
@abstract  Flush buffer to stream.
@param p   Pointer of jfile_t.
@return    0 if success, EOF otherwise.
 */
int jf_flush(jfile_t *p);

int jf_printf(jfile_t *p, const char *fmt, ...);
int jf_putc(int c, jfile_t *p); 
int jf_putc_(int c, jfile_t *p);
int jf_puts(const char *s, jfile_t *p); 
int jf_write(jfile_t *p, char *buf, size_t len);

/*!@func
@abstract  Close jfile_t but does not destroy it.
@param p   Pointer of jfile_t.
@return    0 if success, EOF otherwise.
@note      Even fail, the jfile_t will still be set to not open.
 */
int jf_close(jfile_t *p);

/*!@func
@abstract  Remove file.
@param p   Pointer of jfile_t to be removed.
@return    1 if success, 0 if file does not exist, -1 if failure.
 */
int jf_remove(jfile_t *p);

/*!@func
@abstract  Remove all files in array.
@param fs  Pointer of array of jfile_t to be removed.
@param n   Size of array.
@return    Num of files have been removed if no error, -1 otherwise.
 */
int jf_remove_all(jfile_t **fs, const int n);

/* 
* File Functions
 */

/*!@func
@abstract  Join together two pathes.
@param p1  Pointer to the first path.
@param p2  Pointer to the second path.
@return    Pointer to the joined path if success, NULL otherwise.
@note      Only works for Unix system as the path seperator used in this function is '/'.
 */
char* join_path(const char *p1, const char *p2);

/*!@func
@abstract      Merge several files into one.
@param in_fn   Names of input files to be merged.
@param n       Num of input files.
@param out_fn  Name of output file.
@return        Num of files have been merged.
 */
int merge_files(char **in_fn, const int n, const char *out_fn);

/*!@func
@abstract  Remove file.
@param fn  Name of file to be removed.
@return    1 if success, 0 if file does not exist, -1 if failure.
 */
int remove_file(char *fn);

/*!@func
@abstract  Remove several files.
@param fn  Names of files to be removed.
@param n   Num of files.
@return    Num of files have been removed if no error, -1 otherwise.
 */
int remove_files(char **fn, const int n);

#endif